#!/usr/bin/env python3
"""
Enhanced Plant P450 Enzyme N-terminal/C-terminal Amino Acid Functional Classification Comparative Analysis Script
Features:
1. Extract protein sequence embedding vectors using ESM2
2. Compare classification performance of N-terminal, C-terminal, middle region, random fragment, and full-length sequences
3. Validate signal specificity (avoid false positives from non-relevant sequences)
4. Provide three dimensionality reduction visualizations: PCA, t-SNE, UMAP
5. Export all visualization raw data as CSV files
6. Generate preliminary visualization charts
7. Multi-trial random fragment analysis for robustness validation
8. Add Shuffled N-terminal as additional control
9. Include RandomForest and LogisticRegression classifiers
"""

import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import Counter
import warnings
import random
warnings.filterwarnings('ignore')

# Load ESM2 model using transformers
from transformers import AutoTokenizer, AutoModel

# Machine learning libraries
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.decomposition import PCA
from sklearn.svm import SVC, LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

# Dimensionality reduction methods
try:
    from sklearn.manifold import TSNE
    TSNE_AVAILABLE = True
except ImportError:
    TSNE_AVAILABLE = False
    print("Warning: t-SNE not available, please install scikit-learn")

try:
    import umap.umap_ as umap
    UMAP_AVAILABLE = True
except ImportError:
    UMAP_AVAILABLE = False
    print("Warning: UMAP not available, please install umap-learn")

# Statistical tests
from scipy import stats
import scipy.stats as st

# Progress bar
from tqdm import tqdm

class EnhancedP450Analyzer:
    """Enhanced P450 Analyzer - Includes N-terminal, C-terminal, middle region and random fragment comparative analysis"""
    
    def __init__(self, model_name: str = "facebook/esm2_t12_35M_UR50D", 
                 n_aa: int = 100, 
                 c_aa: int = 100,
                 middle_aa: int = 100,
                 random_aa: int = 100,
                 device: str = "cpu",
                 use_tsne: bool = True,
                 use_umap: bool = True,
                 num_random_trials: int = 20):
        """
        Initialize analyzer
        
        Args:
            model_name: ESM2 model name
            n_aa: Number of N-terminal amino acids to analyze
            c_aa: Number of C-terminal amino acids to analyze
            middle_aa: Number of middle region amino acids to analyze
            random_aa: Number of random fragment amino acids to analyze
            device: Computing device (cpu/cuda)
            use_tsne: Whether to compute t-SNE dimensionality reduction
            use_umap: Whether to compute UMAP dimensionality reduction
            num_random_trials: Number of trials for multi-trial random fragment analysis
        """
        self.n_aa = n_aa
        self.c_aa = c_aa
        self.middle_aa = middle_aa
        self.random_aa = random_aa
        self.device = device if torch.cuda.is_available() and device == "cuda" else "cpu"
        self.use_tsne = use_tsne and TSNE_AVAILABLE
        self.use_umap = use_umap and UMAP_AVAILABLE
        self.num_random_trials = num_random_trials
        self.random_seed = 27  # Fixed random seed for reproducibility
        
        # Load ESM2 model using transformers
        print(f"Loading ESM2 model: {model_name}...")
        try:
            self.tokenizer = AutoTokenizer.from_pretrained(model_name)
            self.model = AutoModel.from_pretrained(model_name)
            self.model = self.model.to(self.device).eval()
            self.embedding_dim = self.model.embeddings.word_embeddings.embedding_dim
            print(f"Model loaded successfully, using device: {self.device}")
            print(f"Embedding dimension: {self.embedding_dim}")
        except Exception as e:
            print(f"Model loading failed: {e}")
            print("Trying smaller model...")
            self.tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
            self.model = AutoModel.from_pretrained("facebook/esm2_t6_8M_UR50D")
            self.model = self.model.to(self.device).eval()
            self.embedding_dim = self.model.embeddings.word_embeddings.embedding_dim
            print(f"Alternative model loaded, embedding dimension: {self.embedding_dim}")
    
    def detect_csv_delimiter(self, filepath: str) -> str:
        """Automatically detect CSV file delimiter"""
        with open(filepath, 'r', encoding='utf-8') as f:
            sample = f.read(4096)
            delimiters = [',', '\t', ';', '|']
            delimiter_counts = {}
            
            for delimiter in delimiters:
                lines = sample.split('\n')[:10]
                counts = [line.count(delimiter) for line in lines if line.strip()]
                if counts:
                    delimiter_counts[delimiter] = max(set(counts), key=counts.count)
            
            if delimiter_counts:
                best_delimiter = max(delimiter_counts, key=delimiter_counts.get)
                print(f"Detected delimiter: '{best_delimiter}'")
                return best_delimiter
            return ','
    
    def load_data(self, fasta_file: str, cluster_file: str) -> Tuple[Dict, Dict]:
        """Load FASTA file and cluster labels"""
        # Load cluster file
        print("Loading cluster label file...")
        delimiter = self.detect_csv_delimiter(cluster_file)
        
        try:
            cluster_df = pd.read_csv(cluster_file, sep=delimiter, engine='python')
        except:
            # Try different delimiters
            for sep in [',', '\t', ';', '|']:
                try:
                    cluster_df = pd.read_csv(cluster_file, sep=sep, engine='python')
                    break
                except:
                    continue
        
        # Parse gene to cluster mapping
        gene_cluster_dict = {}
        for _, row in cluster_df.iterrows():
            cluster_col = 'Cluster' if 'Cluster' in cluster_df.columns else cluster_df.columns[0]
            gene_list_col = 'Gene_List' if 'Gene_List' in cluster_df.columns else cluster_df.columns[2]
            
            cluster_name = str(row[cluster_col]).strip()
            gene_list_str = str(row[gene_list_col]).strip()
            
            # Split gene list
            for sep in [';', ',', '|', '\t', ' ']:
                if sep in gene_list_str:
                    genes = [g.strip() for g in gene_list_str.split(sep) if g.strip()]
                    break
            else:
                genes = [g.strip() for g in gene_list_str.split() if g.strip()]
            
            for gene in genes:
                gene_cluster_dict[gene] = cluster_name
        
        print(f"Loaded cluster labels for {len(gene_cluster_dict)} genes")
        
        # Load FASTA file
        print("\nLoading FASTA file...")
        sequences_dict = {}
        sequence_lengths = []
        
        with open(fasta_file, 'r') as f:
            current_gene = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_gene and current_seq:
                        full_seq = ''.join(current_seq)
                        sequences_dict[current_gene] = full_seq
                        sequence_lengths.append(len(full_seq))
                    
                    # Extract gene name
                    header = line[1:].strip()
                    if '|' in header:
                        current_gene = header.split('|')[0]
                    elif ' ' in header:
                        current_gene = header.split(' ')[0]
                    elif '\t' in header:
                        current_gene = header.split('\t')[0]
                    else:
                        current_gene = header
                    
                    current_seq = []
                else:
                    if line:
                        current_seq.append(line)
            
            if current_gene and current_seq:
                full_seq = ''.join(current_seq)
                sequences_dict[current_gene] = full_seq
                sequence_lengths.append(len(full_seq))
        
        print(f"Loaded {len(sequences_dict)} protein sequences")
        print(f"Sequence length statistics: mean={np.mean(sequence_lengths):.1f}, min={min(sequence_lengths)}, max={max(sequence_lengths)}")
        
        # Align data
        common_genes = set(sequences_dict.keys()) & set(gene_cluster_dict.keys())
        sequences_dict = {g: sequences_dict[g] for g in common_genes}
        gene_cluster_dict = {g: gene_cluster_dict[g] for g in common_genes}
        
        print(f"\nAfter alignment: {len(common_genes)} common genes")
        
        return sequences_dict, gene_cluster_dict
    
    def handle_class_imbalance(self, labels: List) -> Dict:
        """Handle class imbalance"""
        label_counts = Counter(labels)
        print("\n=== Class Distribution Statistics ===")
        for label, count in sorted(label_counts.items(), key=lambda x: x[1], reverse=True):
            print(f"{label}: {count} samples ({count/len(labels)*100:.1f}%)")
        
        # Merge rare classes (samples < 5)
        rare_labels = [label for label, count in label_counts.items() if count < 5]
        if rare_labels:
            print(f"\nRare classes found (samples < 5): {len(rare_labels)} classes")
            labels_processed = ['Rare_Class' if label in rare_labels else label for label in labels]
            new_counts = Counter(labels_processed)
            print("Class distribution after merging rare classes:")
            for label, count in sorted(new_counts.items(), key=lambda x: x[1], reverse=True):
                print(f"{label}: {count} samples ({count/len(labels_processed)*100:.1f}%)")
        else:
            labels_processed = labels
            new_counts = label_counts
        
        # Calculate class balance index
        if len(new_counts) > 1:
            max_count = max(new_counts.values())
            min_count = min(new_counts.values())
            balance_index = min_count / max_count if max_count > 0 else 0
            print(f"Class balance index: {balance_index:.3f} (1.0 indicates perfect balance)")
        
        return {
            'labels': labels_processed,
            'original_labels': labels,
            'rare_labels': rare_labels if rare_labels else [],
            'label_counts': new_counts
        }
    
    def extract_esm2_embeddings(self, sequences: List[str]) -> np.ndarray:
        """Extract ESM2 embedding vectors"""
        print(f"Extracting ESM2 embeddings for {len(sequences)} sequences...")
        all_embeddings = []
        batch_size = 2  # Reduce batch size for large models
        
        with torch.no_grad():
            for i in tqdm(range(0, len(sequences), batch_size), desc="Extracting embeddings"):
                batch_sequences = sequences[i:i+batch_size]
                
                inputs = self.tokenizer(
                    batch_sequences, 
                    padding=True, 
                    truncation=True, 
                    max_length=1024,
                    return_tensors="pt"
                )
                
                inputs = {k: v.to(self.device) for k, v in inputs.items()}
                outputs = self.model(**inputs)
                last_hidden_state = outputs.last_hidden_state
                
                for j in range(len(batch_sequences)):
                    attention_mask = inputs['attention_mask'][j]
                    seq_length = attention_mask.sum().item()
                    
                    seq_embeddings = last_hidden_state[j, 1:seq_length-1]
                    if len(seq_embeddings) > 0:
                        mean_embedding = seq_embeddings.mean(dim=0).cpu().numpy()
                    else:
                        seq_embeddings = last_hidden_state[j, :seq_length]
                        mean_embedding = seq_embeddings.mean(dim=0).cpu().numpy()
                    
                    all_embeddings.append(mean_embedding)
        
        X = np.array(all_embeddings)
        print(f"Feature matrix shape: {X.shape}")
        return X
    
    def shuffle_sequence(self, seq: str) -> str:
        """Shuffle amino acid sequence while preserving amino acid composition"""
        if len(seq) <= 1:
            return seq
        
        # Convert to list for shuffling
        seq_list = list(seq)
        random.shuffle(seq_list)
        return ''.join(seq_list)
    
    def extract_embeddings(self, sequences_dict: Dict, gene_cluster_dict: Dict, 
                          sequence_type: str = 'n_terminal', 
                          random_seed_override: int = None) -> Tuple[np.ndarray, np.ndarray, Dict]:
        """Extract sequence embeddings with safe random fragment generation and shuffled sequences"""
        print(f"\nExtracting ESM2 embeddings for {sequence_type}...")
        
        genes = list(sequences_dict.keys())
        sequences = []
        labels = []
        sequence_info = []
        
        # Track statistics for random fragments
        random_fragment_stats = {
            'full_sequence_used': 0,
            'random_fragments': 0,
            'short_sequences': 0,
            'positions': []  # For position distribution analysis
        }
        
        # Set random seed for reproducibility
        current_random_seed = random_seed_override if random_seed_override is not None else self.random_seed
        np.random.seed(current_random_seed)
        random.seed(current_random_seed)
        
        for gene in genes:
            seq = sequences_dict[gene]
            seq_to_use = None
            
            if sequence_type == 'n_terminal':
                if len(seq) >= self.n_aa:
                    seq_to_use = seq[:self.n_aa]
                else:
                    seq_to_use = seq
            
            elif sequence_type == 'shuffled_n_terminal':
                # Extract N-terminal then shuffle
                if len(seq) >= self.n_aa:
                    n_term_seq = seq[:self.n_aa]
                    seq_to_use = self.shuffle_sequence(n_term_seq)
                else:
                    seq_to_use = self.shuffle_sequence(seq)
            
            elif sequence_type == 'c_terminal':
                if len(seq) >= self.c_aa:
                    seq_to_use = seq[-self.c_aa:]
                else:
                    seq_to_use = seq
            
            elif sequence_type == 'middle':
                if len(seq) >= self.middle_aa:
                    start_idx = (len(seq) - self.middle_aa) // 2
                    seq_to_use = seq[start_idx:start_idx + self.middle_aa]
                else:
                    seq_to_use = seq
            
            elif sequence_type == 'random':
                # Safely generate random fragment
                if len(seq) > self.random_aa:
                    # Normal case: random starting position
                    max_start = len(seq) - self.random_aa
                    start_idx = random.randint(0, max_start)
                    seq_to_use = seq[start_idx:start_idx + self.random_aa]
                    random_fragment_stats['random_fragments'] += 1
                    random_fragment_stats['positions'].append(start_idx / len(seq))
                else:
                    # Sequence too short: use full sequence and mark
                    seq_to_use = seq
                    random_fragment_stats['full_sequence_used'] += 1
                    random_fragment_stats['short_sequences'] += 1
            
            else:  # full_length
                seq_to_use = seq
            
            sequences.append(seq_to_use)
            labels.append(gene_cluster_dict[gene])
            sequence_info.append({
                'gene': gene,
                'length': len(seq),
                'region': sequence_type,
                'region_length': len(seq_to_use),
                'is_full_for_random': (sequence_type == 'random' and len(seq) <= self.random_aa)
            })
        
        # Print random fragment statistics
        if sequence_type == 'random' and random_fragment_stats['positions']:
            print(f"Random fragment statistics:")
            print(f"  - Sequences with random fragments: {random_fragment_stats['random_fragments']}")
            print(f"  - Sequences using full length (too short): {random_fragment_stats['full_sequence_used']}")
            if random_fragment_stats['positions']:
                positions = random_fragment_stats['positions']
                print(f"  - Average normalized position: {np.mean(positions):.3f} (0=N-term, 1=C-term)")
                print(f"  - Position std: {np.std(positions):.3f}")
        
        label_info = self.handle_class_imbalance(labels)
        labels = label_info['labels']
        
        self.label_encoder = LabelEncoder()
        y = self.label_encoder.fit_transform(labels)
        X = self.extract_esm2_embeddings(sequences)
        
        # Add sequence information to label info
        label_info['sequence_info'] = sequence_info
        if sequence_type == 'random':
            label_info['random_fragment_stats'] = random_fragment_stats
        
        return X, y, label_info
    
    def compute_dimension_reduction(self, X: np.ndarray, y: np.ndarray, 
                                   region_name: str) -> Dict:
        """
        Compute 2D coordinates using multiple dimensionality reduction methods
        
        Returns:
            Dictionary containing coordinates from various dimensionality reduction methods
        """
        print(f"\nComputing dimensionality reduction coordinates for {region_name}...")
        
        # Standardization
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        results = {
            'region': region_name,
            'original_features': X_scaled,
            'labels': y,
            'label_names': self.label_encoder.inverse_transform(y)
        }
        
        # 1. PCA dimensionality reduction (for classification)
        pca_full = PCA(n_components=0.95, random_state=42)
        X_pca_full = pca_full.fit_transform(X_scaled)
        results['X_pca'] = X_pca_full
        results['pca_variance_ratio'] = pca_full.explained_variance_ratio_
        
        # 2. PCA 2D (for visualization)
        pca_2d = PCA(n_components=2, random_state=42)
        X_pca_2d = pca_2d.fit_transform(X_scaled)
        results['pca_2d'] = X_pca_2d
        results['pca_explained_var'] = pca_2d.explained_variance_ratio_
        
        print(f"PCA 2D - Explained variance: {pca_2d.explained_variance_ratio_.sum():.3f}")
        
        # 3. t-SNE 2D (if available)
        if self.use_tsne:
            try:
                perplexity = min(30, X_scaled.shape[0] - 1)
                if perplexity >= 5:
                    tsne = TSNE(n_components=2, random_state=42, 
                               perplexity=perplexity, max_iter=1000)
                    X_tsne_2d = tsne.fit_transform(X_scaled)
                    results['tsne_2d'] = X_tsne_2d
                    print("t-SNE computation completed")
                else:
                    print("Insufficient data, skipping t-SNE")
            except Exception as e:
                print(f"t-SNE computation failed: {e}")
        
        # 4. UMAP 2D (if available)
        if self.use_umap:
            try:
                n_neighbors = min(15, X_scaled.shape[0] - 1)
                umap_reducer = umap.UMAP(n_components=2, random_state=42,
                                        n_neighbors=n_neighbors, min_dist=0.1)
                X_umap_2d = umap_reducer.fit_transform(X_scaled)
                results['umap_2d'] = X_umap_2d
                print("UMAP computation completed")
            except Exception as e:
                print(f"UMAP computation failed: {e}")
        
        return results
    
    def compare_feature_regions(self, sequences_dict: Dict, 
                               gene_cluster_dict: Dict) -> Dict:
        """Compare classification performance of different sequence regions"""
        results = {}
        
        # Add middle region, random fragment, and shuffled N-terminal as controls
        regions_to_compare = ['n_terminal', 'shuffled_n_terminal', 'c_terminal', 'middle', 'random', 'full_length']
        region_names = {
            'n_terminal': 'N-terminal',
            'shuffled_n_terminal': 'Shuffled N-terminal',
            'c_terminal': 'C-terminal',
            'middle': 'Middle Region',
            'random': 'Random Fragment',
            'full_length': 'Full-length'
        }
        
        for region in regions_to_compare:
            print(f"\n{'='*60}")
            print(f"Analyzing region: {region_names[region]}")
            
            X, y, label_info = self.extract_embeddings(sequences_dict, gene_cluster_dict, region)
            
            # Compute dimensionality reduction coordinates
            dim_reduction_results = self.compute_dimension_reduction(X, y, region)
            X_pca = dim_reduction_results['X_pca']
            
            print(f"After PCA dimensionality reduction: {X_pca.shape[1]} dimensions")
            print(f"Explained variance ratio: {dim_reduction_results['pca_variance_ratio'].sum():.3f}")
            
            # Evaluate classifiers
            classifiers = {
                'LinearSVM': LinearSVC(C=1.0, class_weight='balanced', random_state=42, max_iter=10000),
                'LogisticRegression': LogisticRegression(penalty='l2', C=0.1, 
                                                       class_weight='balanced', 
                                                       max_iter=10000, random_state=42),
                'RandomForest': RandomForestClassifier(n_estimators=100, max_depth=10,
                                                      min_samples_split=5, min_samples_leaf=2,
                                                      class_weight='balanced', random_state=42)
            }
            
            cv_scores = {}
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            
            for name, clf in classifiers.items():
                scores = cross_val_score(clf, X_pca, y, cv=skf, 
                                        scoring='accuracy', n_jobs=-1)
                cv_scores[name] = {
                    'mean': scores.mean(),
                    'std': scores.std(),
                    'scores': scores
                }
                print(f"{name}: Accuracy = {scores.mean():.4f} ± {scores.std():.4f}")
            
            # Bootstrap validation
            print("Performing Bootstrap validation...")
            bootstrap_scores = []
            n_bootstrap = 100
            n_samples = X_pca.shape[0]
            
            for _ in tqdm(range(n_bootstrap), desc="Bootstrap iterations"):
                indices = np.random.choice(n_samples, n_samples, replace=True)
                X_boot = X_pca[indices]
                y_boot = y[indices]
                
                # Check sample distribution
                unique, counts = np.unique(y_boot, return_counts=True)
                if counts.min() >= 2:
                    try:
                        X_train, X_test, y_train, y_test = train_test_split(
                            X_boot, y_boot, test_size=0.2, stratify=y_boot, random_state=42
                        )
                    except:
                        X_train, X_test, y_train, y_test = train_test_split(
                            X_boot, y_boot, test_size=0.2, random_state=42
                        )
                else:
                    X_train, X_test, y_train, y_test = train_test_split(
                        X_boot, y_boot, test_size=0.2, random_state=42
                    )
                
                best_clf_name = max(cv_scores, key=lambda x: cv_scores[x]['mean'])
                best_clf = classifiers[best_clf_name]
                best_clf.fit(X_train, y_train)
                score = best_clf.score(X_test, y_test)
                bootstrap_scores.append(score)
            
            # Calculate confidence interval
            if len(bootstrap_scores) > 1:
                ci_lower, ci_upper = st.t.interval(0.95, len(bootstrap_scores)-1, 
                                                  loc=np.mean(bootstrap_scores), 
                                                  scale=st.sem(bootstrap_scores))
            else:
                ci_lower = ci_upper = np.mean(bootstrap_scores)
            
            print(f"Bootstrap 95% Confidence Interval: [{ci_lower:.4f}, {ci_upper:.4f}]")
            
            # Store results
            results[region] = {
                'X': X_pca,
                'y': y,
                'cv_scores': cv_scores,
                'bootstrap_scores': bootstrap_scores,
                'bootstrap_ci': (ci_lower, ci_upper),
                'label_info': label_info,
                'dim_reduction': dim_reduction_results,
                'best_classifier': max(cv_scores, key=lambda x: cv_scores[x]['mean']),
                'region_name': region_names[region]
            }
        
        return results
    
    def evaluate_multiple_random_trials(self, sequences_dict: Dict, gene_cluster_dict: Dict, 
                                       n_trials: int = None) -> Dict:
        """
        Perform multiple random fragment trials for robustness validation
        
        Args:
            sequences_dict: Dictionary of gene to sequence mappings
            gene_cluster_dict: Dictionary of gene to cluster mappings
            n_trials: Number of trials (defaults to self.num_random_trials)
        
        Returns:
            Dictionary containing trial results and statistics
        """
        if n_trials is None:
            n_trials = self.num_random_trials
            
        print(f"\n{'='*60}")
        print(f"Performing {n_trials} random fragment trials for robustness validation")
        print(f"{'='*60}")
        
        trial_results = {
            'accuracies': [],
            'positions_all': [],  # Track position distribution across all trials
            'trial_details': []
        }
        
        for trial in range(n_trials):
            print(f"\nTrial {trial+1}/{n_trials}")
            
            # Use different random seed for each trial
            random_seed = 42 + trial * 1000  # Different seed for each trial
            
            try:
                # Extract current trial's random fragment embeddings
                X_random, y_random, label_info = self.extract_embeddings(
                    sequences_dict, gene_cluster_dict, 'random', random_seed_override=random_seed
                )
                
                # Compute dimensionality reduction
                dim_results = self.compute_dimension_reduction(X_random, y_random, f'random_trial_{trial}')
                X_pca = dim_results['X_pca']
                
                # Evaluate with linear SVM (consistent with main analysis)
                clf = LinearSVC(C=1.0, class_weight='balanced', random_state=42, max_iter=10000)
                scores = cross_val_score(clf, X_pca, y_random, cv=5, 
                                       scoring='accuracy', n_jobs=-1)
                
                trial_acc = scores.mean()
                trial_results['accuracies'].append(trial_acc)
                
                # Record trial details
                trial_detail = {
                    'trial': trial + 1,
                    'accuracy': trial_acc,
                    'accuracy_scores': scores.tolist(),
                    'random_seed': random_seed
                }
                
                # Add random fragment statistics if available
                if 'random_fragment_stats' in label_info:
                    stats = label_info['random_fragment_stats']
                    trial_detail['stats'] = stats
                    trial_results['positions_all'].extend(stats['positions'])
                
                trial_results['trial_details'].append(trial_detail)
                
                print(f"  Trial accuracy: {trial_acc:.4f}")
                
            except Exception as e:
                print(f"  Trial {trial+1} failed: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        # Calculate overall statistics
        if trial_results['accuracies']:
            accuracies = trial_results['accuracies']
            mean_acc = np.mean(accuracies)
            std_acc = np.std(accuracies)
            
            # Calculate confidence interval
            if len(accuracies) > 1:
                ci_lower, ci_upper = st.t.interval(0.95, len(accuracies)-1,
                                                  loc=mean_acc, scale=st.sem(accuracies))
            else:
                ci_lower = ci_upper = mean_acc
            
            print(f"\n{'='*60}")
            print("MULTI-TRIAL RANDOM FRAGMENT ANALYSIS RESULTS:")
            print(f"{'='*60}")
            print(f"Number of successful trials: {len(accuracies)}")
            print(f"Mean accuracy: {mean_acc:.4f}")
            print(f"Standard deviation: {std_acc:.4f}")
            print(f"95% Confidence Interval: [{ci_lower:.4f}, {ci_upper:.4f}]")
            print(f"Accuracy range: [{min(accuracies):.4f}, {max(accuracies):.4f}]")
            
            # Position distribution analysis
            if trial_results['positions_all']:
                positions = trial_results['positions_all']
                from scipy.stats import kstest
                stat, p_value = kstest(positions, 'uniform')
                print(f"\nRandom Position Uniformity Test:")
                print(f"  K-S test statistic: {stat:.4f}")
                print(f"  p-value: {p_value:.4f}")
                if p_value > 0.05:
                    print("  ✓ Uniform distribution cannot be rejected")
                else:
                    print("  ⚠️ Significantly non-uniform distribution")
        
        return trial_results
    
    def plot_random_trial_distributions(self, trial_results: Dict, output_dir: str = "enhanced_p450_region_analysis"):
        """Plot multi-trial random fragment analysis distributions"""
        # Create plot directory
        plot_dir = Path(output_dir) / "plots"
        plot_dir.mkdir(exist_ok=True, parents=True)
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # 1. Accuracy distribution
        accuracies = trial_results['accuracies']
        axes[0].hist(accuracies, bins=10, alpha=0.7, color='steelblue', edgecolor='black')
        axes[0].axvline(np.mean(accuracies), color='red', linestyle='--', 
                        label=f'Mean: {np.mean(accuracies):.4f}')
        axes[0].set_xlabel('Accuracy')
        axes[0].set_ylabel('Frequency')
        axes[0].set_title(f'Random Fragment Accuracy Distribution (n={len(accuracies)} trials)')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # 2. Position distribution (if available)
        if trial_results['positions_all']:
            positions = trial_results['positions_all']
            axes[1].hist(positions, bins=20, alpha=0.7, color='forestgreen', 
                        edgecolor='black', density=True)
            axes[1].axhline(y=1.0, color='red', linestyle='--', 
                           label='Uniform Distribution')
            axes[1].set_xlabel('Normalized Start Position (0=N-term, 1=C-term)')
            axes[1].set_ylabel('Density')
            axes[1].set_title('Distribution of Random Fragment Start Positions')
            axes[1].legend()
            axes[1].grid(True, alpha=0.3)
        else:
            axes[1].text(0.5, 0.5, 'Position data not available', 
                        ha='center', va='center', transform=axes[1].transAxes)
            axes[1].set_title('Position Distribution')
        
        plt.tight_layout()
        plt.savefig(f"{plot_dir}/random_fragment_distributions.png", 
                    dpi=300, bbox_inches='tight')
        plt.savefig(f"{plot_dir}/random_fragment_distributions.pdf", 
                    bbox_inches='tight')
        plt.show()
    
    def export_visualization_data(self, results: Dict, output_dir: str):
        """Export all visualization raw data"""
        print(f"\n{'='*60}")
        print("Exporting visualization raw data...")
        Path(output_dir).mkdir(exist_ok=True, parents=True)
        
        # 1. Export dimensionality reduction coordinates
        for region, region_data in results.items():
            dim_data = region_data['dim_reduction']
            y = region_data['y']
            label_names = dim_data['label_names']
            
            # Create dataframe
            plot_df = pd.DataFrame({
                'Sample_ID': range(len(y)),
                'Numeric_Label': y,
                'Class_Label': label_names
            })
            
            # Add PCA coordinates
            if 'pca_2d' in dim_data:
                pca_2d = dim_data['pca_2d']
                plot_df['PCA_Dim1'] = pca_2d[:, 0]
                plot_df['PCA_Dim2'] = pca_2d[:, 1]
                plot_df['PCA_Var1'] = dim_data['pca_explained_var'][0]
                plot_df['PCA_Var2'] = dim_data['pca_explained_var'][1]
            
            # Add t-SNE coordinates
            if 'tsne_2d' in dim_data:
                tsne_2d = dim_data['tsne_2d']
                plot_df['tSNE_Dim1'] = tsne_2d[:, 0]
                plot_df['tSNE_Dim2'] = tsne_2d[:, 1]
            
            # Add UMAP coordinates
            if 'umap_2d' in dim_data:
                umap_2d = dim_data['umap_2d']
                plot_df['UMAP_Dim1'] = umap_2d[:, 0]
                plot_df['UMAP_Dim2'] = umap_2d[:, 1]
            
            # Save to CSV
            region_name = region_data['region_name'].replace('-', '_').replace(' ', '_')
            filename = f"{output_dir}/{region_name}_projection_data.csv"
            plot_df.to_csv(filename, index=False, encoding='utf-8')
            print(f"Saved: {filename}")
        
        # 2. Export classifier performance data
        performance_data = []
        for region, region_data in results.items():
            region_name = region_data['region_name']
            for clf_name, clf_data in region_data['cv_scores'].items():
                performance_data.append({
                    'Region': region_name,
                    'Classifier': clf_name,
                    'Accuracy_Mean': clf_data['mean'],
                    'Accuracy_SD': clf_data['std'],
                    'Accuracy_Values': ';'.join([f"{x:.4f}" for x in clf_data['scores']])
                })
        
        perf_df = pd.DataFrame(performance_data)
        perf_df.to_csv(f"{output_dir}/classifier_performance.csv", index=False, encoding='utf-8')
        print(f"Saved: {output_dir}/classifier_performance.csv")
        
        # 3. Export detailed accuracy data for R visualization
        detailed_performance_data = []
        for region, region_data in results.items():
            region_name = region_data['region_name']
            for clf_name, clf_data in region_data['cv_scores'].items():
                for fold_idx, score in enumerate(clf_data['scores']):
                    detailed_performance_data.append({
                        'Region': region_name,
                        'Classifier': clf_name,
                        'Fold': fold_idx + 1,
                        'Accuracy': score
                    })
        
        detailed_df = pd.DataFrame(detailed_performance_data)
        detailed_df.to_csv(f"{output_dir}/detailed_accuracy_data.csv", index=False, encoding='utf-8')
        print(f"Saved: {output_dir}/detailed_accuracy_data.csv")
        
        # 4. Export Bootstrap data
        bootstrap_data = []
        for region, region_data in results.items():
            region_name = region_data['region_name']
            for i, score in enumerate(region_data['bootstrap_scores']):
                bootstrap_data.append({
                    'Region': region_name,
                    'Iteration': i,
                    'Accuracy': score
                })
        
        bootstrap_df = pd.DataFrame(bootstrap_data)
        bootstrap_df.to_csv(f"{output_dir}/bootstrap_distribution.csv", index=False, encoding='utf-8')
        print(f"Saved: {output_dir}/bootstrap_distribution.csv")
        
        # 5. Export class distribution data
        if 'n_terminal' in results:
            label_counts = results['n_terminal']['label_info']['label_counts']
            label_df = pd.DataFrame(list(label_counts.items()), 
                                   columns=['Class', 'Count'])
            label_df.to_csv(f"{output_dir}/class_distribution.csv", index=False, encoding='utf-8')
            print(f"Saved: {output_dir}/class_distribution.csv")
        
        # 6. Export inter-region performance comparison data
        comparison_data = []
        if 'n_terminal' in results and 'shuffled_n_terminal' in results:
            n_term_scores = results['n_terminal']['cv_scores']['LinearSVM']['scores']
            shuffled_scores = results['shuffled_n_terminal']['cv_scores']['LinearSVM']['scores']
            
            for i, (n_score, s_score) in enumerate(zip(n_term_scores, shuffled_scores)):
                comparison_data.append({
                    'Fold': i+1,
                    'N_terminal_Accuracy': n_score,
                    'Shuffled_N_terminal_Accuracy': s_score,
                    'Difference_N_Shuffled': n_score - s_score
                })
        
        # Add other region comparisons
        if 'c_terminal' in results and 'middle' in results and 'random' in results:
            c_term_scores = results['c_terminal']['cv_scores']['LinearSVM']['scores']
            middle_scores = results['middle']['cv_scores']['LinearSVM']['scores']
            random_scores = results['random']['cv_scores']['LinearSVM']['scores']
            
            for i, (c_score, m_score, r_score) in enumerate(zip(c_term_scores, middle_scores, random_scores)):
                if i < len(comparison_data):
                    comparison_data[i]['C_terminal_Accuracy'] = c_score
                    comparison_data[i]['Middle_Accuracy'] = m_score
                    comparison_data[i]['Random_Accuracy'] = r_score
        
        if comparison_data:
            comparison_df = pd.DataFrame(comparison_data)
            comparison_df.to_csv(f"{output_dir}/region_performance_comparison.csv", index=False, encoding='utf-8')
            print(f"Saved: {output_dir}/region_performance_comparison.csv")
        
        # 7. Export summary statistics for bar plot in R
        summary_data = []
        for region, region_data in results.items():
            region_name = region_data['region_name']
            for clf_name, clf_data in region_data['cv_scores'].items():
                summary_data.append({
                    'Region': region_name,
                    'Classifier': clf_name,
                    'Mean_Accuracy': clf_data['mean'],
                    'SD_Accuracy': clf_data['std'],
                    'Min_Accuracy': min(clf_data['scores']),
                    'Max_Accuracy': max(clf_data['scores']),
                    'N_Folds': len(clf_data['scores'])
                })
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(f"{output_dir}/summary_statistics.csv", index=False, encoding='utf-8')
        print(f"Saved: {output_dir}/summary_statistics.csv")
        
        print(f"\nAll raw data exported to: {output_dir}/")
        print("File list:")
        for region, region_data in results.items():
            region_name = region_data['region_name'].replace('-', '_').replace(' ', '_')
            print(f"1. {region_name}_projection_data.csv - {region_data['region_name']} dimensionality reduction coordinates (PCA/t-SNE/UMAP)")
        print("4. classifier_performance.csv - Classifier performance comparison")
        print("5. detailed_accuracy_data.csv - Detailed accuracy data for each fold")
        print("6. bootstrap_distribution.csv - Bootstrap accuracy distribution")
        print("7. class_distribution.csv - Class distribution statistics")
        print("8. region_performance_comparison.csv - Detailed inter-region performance comparison")
        print("9. summary_statistics.csv - Summary statistics for bar plots")
    
    def visualize_results(self, results: Dict, output_dir: str = "results"):
        """Generate preliminary visualization charts"""
        print(f"\n{'='*60}")
        print("Generating preliminary visualization charts...")
        
        # First export data
        self.export_visualization_data(results, output_dir)
        
        # Create plot directory
        plot_dir = Path(output_dir) / "plots"
        plot_dir.mkdir(exist_ok=True, parents=True)
        
        # Set plotting style
        plt.style.use('seaborn-v0_8-whitegrid')
        sns.set_palette("husl")
        
        # 1. Classifier performance comparison plot (all regions) - Bar chart
        fig1, ax1 = plt.subplots(figsize=(16, 8))
        regions = []
        classifiers = []
        means = []
        errors = []
        
        # Define colors
        region_colors = {
            'N-terminal': '#2E86AB',
            'Shuffled N-terminal': '#4A6FA5',
            'C-terminal': '#A23B72',
            'Middle Region': '#F18F01',
            'Random Fragment': '#73B761',
            'Full-length': '#C73E1D'
        }
        
        classifier_colors = {
            'LinearSVM': '#FF6B6B',
            'LogisticRegression': '#4ECDC4',
            'RandomForest': '#45B7D1'
        }
        
        # Collect data
        for region, region_data in results.items():
            region_name = region_data['region_name']
            for clf_name, clf_data in region_data['cv_scores'].items():
                regions.append(region_name)
                classifiers.append(clf_name)
                means.append(clf_data['mean'])
                errors.append(clf_data['std'])
        
        # Create grouped bar chart
        n_classifiers = len(set(classifiers))
        n_regions = len(set(regions))
        
        x = np.arange(n_regions)
        width = 0.25
        
        # Group by classifier
        classifier_list = sorted(set(classifiers))
        region_list = sorted(set(regions), key=lambda x: list(region_colors.keys()).index(x) if x in region_colors else 999)
        
        for i, clf in enumerate(classifier_list):
            clf_means = [means[j] for j in range(len(means)) if classifiers[j] == clf]
            clf_errors = [errors[j] for j in range(len(errors)) if classifiers[j] == clf]
            
            # Ensure correct ordering
            ordered_means = []
            ordered_errors = []
            for region in region_list:
                for j, r in enumerate(regions):
                    if r == region and classifiers[j] == clf:
                        ordered_means.append(means[j])
                        ordered_errors.append(errors[j])
                        break
            
            positions = x + i * width - (n_classifiers * width / 2) + width/2
            ax1.bar(positions, ordered_means, width, 
                   label=clf, yerr=ordered_errors, capsize=5, 
                   alpha=0.8, color=classifier_colors.get(clf, '#999999'), 
                   edgecolor='black')
        
        ax1.set_xlabel('Sequence Region')
        ax1.set_ylabel('Accuracy')
        ax1.set_title('Classifier Performance Comparison Across Different Sequence Regions', fontsize=14, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels(region_list, rotation=45, ha='right')
        ax1.legend(title='Classifier', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.set_ylim(0, 1)
        ax1.grid(True, alpha=0.3, axis='y')
        
        # Add reference line: random guessing level
        if 'n_terminal' in results:
            n_classes = len(results['n_terminal']['label_info']['label_counts'])
            random_guess = 1.0 / n_classes
            ax1.axhline(y=random_guess, color='red', linestyle='--', alpha=0.5, 
                       label=f'Random guess ({random_guess:.2f})')
            ax1.legend(title='Classifier', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        plt.savefig(f"{plot_dir}/classifier_performance_comparison_bar.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{plot_dir}/classifier_performance_comparison_bar.pdf", bbox_inches='tight')
        
        # 2. Region performance comparison (averaged across classifiers) - Bar chart
        fig2, ax2 = plt.subplots(figsize=(12, 7))
        
        region_means = {}
        region_stds = {}
        
        for region, region_data in results.items():
            region_name = region_data['region_name']
            region_accuracies = []
            
            for clf_name, clf_data in region_data['cv_scores'].items():
                region_accuracies.extend(clf_data['scores'])
            
            region_means[region_name] = np.mean(region_accuracies)
            region_stds[region_name] = np.std(region_accuracies)
        
        # Sort regions by mean accuracy
        sorted_regions = sorted(region_means.items(), key=lambda x: x[1], reverse=True)
        sorted_region_names = [r[0] for r in sorted_regions]
        sorted_means = [r[1] for r in sorted_regions]
        sorted_stds = [region_stds[r[0]] for r in sorted_regions]
        
        # Create bar positions
        x_pos = np.arange(len(sorted_region_names))
        
        # Plot bars
        bars = ax2.bar(x_pos, sorted_means, yerr=sorted_stds, capsize=5, 
                      alpha=0.8, color=[region_colors.get(r, '#999999') for r in sorted_region_names],
                      edgecolor='black')
        
        # Add value labels on top of bars
        for i, (bar, mean_val) in enumerate(zip(bars, sorted_means)):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{mean_val:.3f}', ha='center', va='bottom', fontsize=10)
        
        ax2.set_xlabel('Sequence Region')
        ax2.set_ylabel('Average Accuracy (across classifiers)')
        ax2.set_title('Sequence Region Classification Performance Comparison', fontsize=14, fontweight='bold')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(sorted_region_names, rotation=45, ha='right')
        ax2.set_ylim(0, 1)
        ax2.grid(True, alpha=0.3, axis='y')
        
        # Highlight N-terminal region
        if 'N-terminal' in sorted_region_names:
            idx = sorted_region_names.index('N-terminal')
            bars[idx].set_hatch('//')
            bars[idx].set_edgecolor('red')
            bars[idx].set_linewidth(2)
        
        plt.tight_layout()
        plt.savefig(f"{plot_dir}/region_performance_comparison_bar.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{plot_dir}/region_performance_comparison_bar.pdf", bbox_inches='tight')
        
        # 3. Heatmap of classifier performance across regions
        fig3, ax3 = plt.subplots(figsize=(12, 8))
        
        # Create performance matrix
        region_list_heatmap = sorted(set(regions), key=lambda x: list(region_colors.keys()).index(x) if x in region_colors else 999)
        classifier_list_heatmap = sorted(set(classifiers))
        
        performance_matrix = np.zeros((len(region_list_heatmap), len(classifier_list_heatmap)))
        
        for i, region in enumerate(region_list_heatmap):
            for j, clf in enumerate(classifier_list_heatmap):
                # Find the corresponding mean accuracy
                for k in range(len(regions)):
                    if regions[k] == region and classifiers[k] == clf:
                        performance_matrix[i, j] = means[k]
                        break
        
        # Create heatmap
        im = ax3.imshow(performance_matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)
        
        # Add text annotations
        for i in range(len(region_list_heatmap)):
            for j in range(len(classifier_list_heatmap)):
                text = ax3.text(j, i, f'{performance_matrix[i, j]:.3f}',
                              ha="center", va="center", color="black", fontsize=10)
        
        ax3.set_xticks(np.arange(len(classifier_list_heatmap)))
        ax3.set_yticks(np.arange(len(region_list_heatmap)))
        ax3.set_xticklabels(classifier_list_heatmap)
        ax3.set_yticklabels(region_list_heatmap)
        ax3.set_xlabel('Classifier')
        ax3.set_ylabel('Sequence Region')
        ax3.set_title('Classification Accuracy Heatmap', fontsize=14, fontweight='bold')
        
        # Add color bar
        cbar = ax3.figure.colorbar(im, ax=ax3)
        cbar.ax.set_ylabel('Accuracy', rotation=-90, va="bottom")
        
        plt.tight_layout()
        plt.savefig(f"{plot_dir}/performance_heatmap.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{plot_dir}/performance_heatmap.pdf", bbox_inches='tight')
        
        # 4. N-terminal vs Shuffled N-terminal comparison
        if 'n_terminal' in results and 'shuffled_n_terminal' in results:
            fig4, axes4 = plt.subplots(1, 2, figsize=(14, 6))
            
            # Left: Bar chart comparison
            n_term_acc = results['n_terminal']['cv_scores']['LinearSVM']['mean']
            shuffled_acc = results['shuffled_n_terminal']['cv_scores']['LinearSVM']['mean']
            n_term_std = results['n_terminal']['cv_scores']['LinearSVM']['std']
            shuffled_std = results['shuffled_n_terminal']['cv_scores']['LinearSVM']['std']
            
            bars = axes4[0].bar(['N-terminal', 'Shuffled N-terminal'], 
                               [n_term_acc, shuffled_acc],
                               yerr=[n_term_std, shuffled_std],
                               capsize=10,
                               color=['#2E86AB', '#4A6FA5'],
                               alpha=0.8,
                               edgecolor='black')
            
            axes4[0].set_ylabel('Accuracy')
            axes4[0].set_title('N-terminal vs Shuffled N-terminal Performance (LinearSVM)', fontsize=12, fontweight='bold')
            axes4[0].set_ylim(0, 1)
            axes4[0].grid(True, alpha=0.3, axis='y')
            
            # Add value labels
            for i, (bar, acc) in enumerate(zip(bars, [n_term_acc, shuffled_acc])):
                axes4[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                            f'{acc:.3f}', ha='center', va='bottom', fontsize=10)
            
            # Right: Fold-by-fold comparison
            n_term_scores = results['n_terminal']['cv_scores']['LinearSVM']['scores']
            shuffled_scores = results['shuffled_n_terminal']['cv_scores']['LinearSVM']['scores']
            
            x_fold = np.arange(len(n_term_scores))
            width_fold = 0.35
            
            axes4[1].bar(x_fold - width_fold/2, n_term_scores, width_fold, 
                        label='N-terminal', color='#2E86AB', alpha=0.8, edgecolor='black')
            axes4[1].bar(x_fold + width_fold/2, shuffled_scores, width_fold, 
                        label='Shuffled N-terminal', color='#4A6FA5', alpha=0.8, edgecolor='black')
            
            axes4[1].set_xlabel('Fold')
            axes4[1].set_ylabel('Accuracy')
            axes4[1].set_title('Fold-by-Fold Comparison', fontsize=12, fontweight='bold')
            axes4[1].set_xticks(x_fold)
            axes4[1].set_xticklabels([f'Fold {i+1}' for i in range(len(n_term_scores))])
            axes4[1].set_ylim(0, 1)
            axes4[1].legend()
            axes4[1].grid(True, alpha=0.3, axis='y')
            
            plt.tight_layout()
            plt.savefig(f"{plot_dir}/n_terminal_vs_shuffled_comparison.png", dpi=300, bbox_inches='tight')
            plt.savefig(f"{plot_dir}/n_terminal_vs_shuffled_comparison.pdf", bbox_inches='tight')
        
        # 5. Dimensionality reduction visualization (N-terminal, Shuffled N-terminal, C-terminal, middle and random fragments)
        regions_to_plot = ['n_terminal', 'shuffled_n_terminal', 'c_terminal', 'middle', 'random']
        region_names_plot = {
            'n_terminal': 'N-terminal',
            'shuffled_n_terminal': 'Shuffled N-terminal',
            'c_terminal': 'C-terminal',
            'middle': 'Middle Region',
            'random': 'Random Fragment'
        }
        
        if all(r in results for r in regions_to_plot):
            fig5, axes = plt.subplots(len(regions_to_plot), 3, figsize=(18, 6*len(regions_to_plot)))
            
            for row_idx, region in enumerate(regions_to_plot):
                dim_data = results[region]['dim_reduction']
                y = results[region]['y']
                label_names = dim_data['label_names']
                
                # PCA plot
                if 'pca_2d' in dim_data:
                    pca_2d = dim_data['pca_2d']
                    var_exp = dim_data['pca_explained_var']
                    
                    scatter = axes[row_idx, 0].scatter(pca_2d[:, 0], pca_2d[:, 1], 
                                                      c=y, cmap='tab20', alpha=0.7, s=50)
                    axes[row_idx, 0].set_xlabel(f'PC1 ({var_exp[0]:.1%})')
                    axes[row_idx, 0].set_ylabel(f'PC2 ({var_exp[1]:.1%})')
                    axes[row_idx, 0].set_title(f'PCA - {region_names_plot[region]}')
                    axes[row_idx, 0].grid(True, alpha=0.3)
                
                # t-SNE plot
                if 'tsne_2d' in dim_data:
                    tsne_2d = dim_data['tsne_2d']
                    axes[row_idx, 1].scatter(tsne_2d[:, 0], tsne_2d[:, 1], 
                                            c=y, cmap='tab20', alpha=0.7, s=50)
                    axes[row_idx, 1].set_xlabel('t-SNE Dimension 1')
                    axes[row_idx, 1].set_ylabel('t-SNE Dimension 2')
                    axes[row_idx, 1].set_title(f't-SNE - {region_names_plot[region]}')
                    axes[row_idx, 1].grid(True, alpha=0.3)
                
                # UMAP plot
                if 'umap_2d' in dim_data:
                    umap_2d = dim_data['umap_2d']
                    scatter = axes[row_idx, 2].scatter(umap_2d[:, 0], umap_2d[:, 1], 
                                                      c=y, cmap='tab20', alpha=0.7, s=50)
                    axes[row_idx, 2].set_xlabel('UMAP Dimension 1')
                    axes[row_idx, 2].set_ylabel('UMAP Dimension 2')
                    axes[row_idx, 2].set_title(f'UMAP - {region_names_plot[region]}')
                    axes[row_idx, 2].grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(f"{plot_dir}/region_dimension_reduction.png", dpi=300, bbox_inches='tight')
            plt.savefig(f"{plot_dir}/region_dimension_reduction.pdf", bbox_inches='tight')
        
        # Show plots
        plt.show()
        print(f"\nPreliminary charts saved to: {plot_dir}/")
        print("For more sophisticated custom charts, use the exported CSV data with R/ggplot2")
        
        # Print summary table
        print(f"\n{'='*80}")
        print("SUMMARY OF CLASSIFICATION ACCURACY")
        print(f"{'='*80}")
        print(f"{'Region':<25} {'LinearSVM':<12} {'LogisticRegression':<20} {'RandomForest':<15}")
        print(f"{'-'*80}")
        
        for region_name in sorted(set(regions), key=lambda x: list(region_colors.keys()).index(x) if x in region_colors else 999):
            svm_acc = ""
            lr_acc = ""
            rf_acc = ""
            
            for i in range(len(regions)):
                if regions[i] == region_name:
                    if classifiers[i] == 'LinearSVM':
                        svm_acc = f"{means[i]:.4f} ± {errors[i]:.4f}"
                    elif classifiers[i] == 'LogisticRegression':
                        lr_acc = f"{means[i]:.4f} ± {errors[i]:.4f}"
                    elif classifiers[i] == 'RandomForest':
                        rf_acc = f"{means[i]:.4f} ± {errors[i]:.4f}"
            
            print(f"{region_name:<25} {svm_acc:<12} {lr_acc:<20} {rf_acc:<15}")

def main():
    """Main function"""
    print("="*80)
    print("Enhanced Plant P450 Enzyme N-terminal/C-terminal/Middle/Random Fragment Amino Acid Functional Classification Comparative Analysis")
    print("Purpose: Validate N-terminal signal specificity, avoid false positives")
    print("="*80)
    
    # Parameter settings
    FASTA_FILE = "P450_unique_pep_final.fasta"  # Replace with your FASTA file path
    CLUSTER_FILE = "P450_unique_pep_final_esm_clusters_genes.csv"  # Replace with your cluster label file path
    OUTPUT_DIR = "enhanced_p450_region_analysis"
    
    # Initialize analyzer
    analyzer = EnhancedP450Analyzer(
        model_name="facebook/esm2_t33_650M_UR50D",  # Use larger model
        n_aa=100,      # Analyze 100 N-terminal amino acids
        c_aa=100,      # Analyze 100 C-terminal amino acids
        middle_aa=100, # Analyze 100 middle region amino acids
        random_aa=100, # Analyze 100 random fragment amino acids
        device="cpu",  # Use CPU (change to "cuda" if GPU available)
        use_tsne=True, # Enable t-SNE
        use_umap=True, # Enable UMAP
        num_random_trials=20  # Number of trials for random fragment analysis
    )
    
    # Load data
    try:
        sequences_dict, gene_cluster_dict = analyzer.load_data(FASTA_FILE, CLUSTER_FILE)
    except Exception as e:
        print(f"Data loading failed: {e}")
        return
    
    if len(sequences_dict) == 0:
        print("Error: No gene sequences found")
        return
    
    # Compare features from different regions
    print("\n" + "="*80)
    print("Starting comparison of classification performance for N-terminal, Shuffled N-terminal, C-terminal, middle region, random fragment, and full-length sequences...")
    print("Shuffled N-terminal, C-terminal, middle region, and random fragment serve as controls to validate N-terminal signal specificity")
    print("="*80)
    
    try:
        results = analyzer.compare_feature_regions(sequences_dict, gene_cluster_dict)
        
        # Statistical tests - multiple comparisons
        print(f"\n{'='*60}")
        print("Statistical Significance Testing:")
        
        # Prepare SVM scores for all regions for comparison
        region_scores = {}
        for region, region_data in results.items():
            region_name = region_data['region_name']
            region_scores[region_name] = region_data['cv_scores']['LinearSVM']['scores']
        
        # Perform paired t-tests for N-terminal vs. other regions
        if 'N-terminal' in region_scores:
            n_term_scores = region_scores['N-terminal']
            
            for other_region, other_scores in region_scores.items():
                if other_region != 'N-terminal' and len(n_term_scores) == len(other_scores):
                    t_stat, p_value = stats.ttest_rel(n_term_scores, other_scores)
                    mean_diff = np.mean(n_term_scores) - np.mean(other_scores)
                    pooled_std = np.sqrt((np.std(n_term_scores, ddof=1)**2 + 
                                        np.std(other_scores, ddof=1)**2) / 2)
                    cohens_d = mean_diff / pooled_std if pooled_std != 0 else 0
                    
                    print(f"\nN-terminal vs {other_region}:")
                    print(f"   Paired t-test: t = {t_stat:.4f}, p = {p_value:.4f}")
                    print(f"   Mean difference: {mean_diff:.4f}")
                    print(f"   Effect size (Cohen's d): {cohens_d:.4f}")
                    
                    # Effect size interpretation
                    if abs(cohens_d) > 0.8:
                        effect_size = "Large effect"
                    elif abs(cohens_d) > 0.5:
                        effect_size = "Medium effect"
                    elif abs(cohens_d) > 0.2:
                        effect_size = "Small effect"
                    else:
                        effect_size = "Negligible effect"
                    
                    print(f"   Effect size interpretation: {effect_size}")
            
            # Multiple comparison correction (Bonferroni)
            n_comparisons = len(region_scores) - 1  # N-terminal vs. all other regions
            alpha = 0.05
            bonferroni_alpha = alpha / n_comparisons
            
            print(f"\nMultiple comparison correction (Bonferroni, α={alpha}):")
            print(f"   Number of comparisons: {n_comparisons}")
            print(f"   Corrected significance level: {bonferroni_alpha:.4f}")
        
        # Multi-trial random fragment analysis for robustness
        print(f"\n{'='*80}")
        print("Performing multi-trial random fragment analysis for robustness validation...")
        print(f"{'='*80}")
        
        random_trial_results = analyzer.evaluate_multiple_random_trials(
            sequences_dict, gene_cluster_dict
        )
        
        # Plot random trial distributions
        analyzer.plot_random_trial_distributions(random_trial_results, OUTPUT_DIR)
        
        # Compare single trial vs multi-trial results
        if 'random' in results and random_trial_results['accuracies']:
            single_trial_acc = results['random']['cv_scores']['LinearSVM']['mean']
            multi_trial_mean = np.mean(random_trial_results['accuracies'])
            multi_trial_std = np.std(random_trial_results['accuracies'])
            
            print(f"\n{'='*60}")
            print("COMPARISON: Single Trial vs Multi-Trial Random Fragment Analysis")
            print(f"{'='*60}")
            print(f"Single trial accuracy: {single_trial_acc:.4f}")
            print(f"Multi-trial mean accuracy: {multi_trial_mean:.4f} ± {multi_trial_std:.4f}")
            print(f"Difference: {single_trial_acc - multi_trial_mean:+.4f}")
            
            # Check if single trial is within 2 standard deviations of multi-trial distribution
            ci_lower = multi_trial_mean - 2 * multi_trial_std
            ci_upper = multi_trial_mean + 2 * multi_trial_std
            
            if ci_lower <= single_trial_acc <= ci_upper:
                print(f"✓ Single trial result is within 2σ of multi-trial distribution.")
                print("  Single trial result is representative of random fragment performance.")
            else:
                print(f"⚠️ Single trial result is outside 2σ of multi-trial distribution.")
                print("  Consider using multi-trial mean for more robust comparison.")
            
            # Update random fragment results with multi-trial statistics if needed
            if abs(single_trial_acc - multi_trial_mean) > 0.05:  # If difference > 5%
                print(f"\nUpdating random fragment results with multi-trial statistics...")
                # Create a synthetic "random_multi_trial" region in results
                results['random_multi_trial'] = {
                    'region_name': 'Random Fragment (Multi-trial)',
                    'cv_scores': {
                        'LinearSVM': {
                            'mean': multi_trial_mean,
                            'std': multi_trial_std,
                            'scores': random_trial_results['accuracies']
                        }
                    },
                    'multi_trial_stats': {
                        'mean': multi_trial_mean,
                        'std': multi_trial_std,
                        'ci': (ci_lower, ci_upper),
                        'n_trials': len(random_trial_results['accuracies'])
                    }
                }
        
        # Visualize and export data
        analyzer.visualize_results(results, OUTPUT_DIR)
        
        # Print key conclusions
        print(f"\n{'='*80}")
        print("Key Conclusions:")
        
        # Extract best classifier accuracy for each region
        region_accuracies = {}
        for region, region_data in results.items():
            region_name = region_data['region_name']
            best_clf_name = region_data['best_classifier']
            accuracy = region_data['cv_scores'][best_clf_name]['mean']
            region_accuracies[region_name] = accuracy
            print(f"{region_name} classification accuracy: {accuracy:.4f} ({best_clf_name})")
        
        # Find region with best performance
        best_region = max(region_accuracies, key=region_accuracies.get)
        print(f"\nBest performing region: {best_region} (accuracy: {region_accuracies[best_region]:.4f})")
        
        # Signal specificity analysis
        print(f"\n{'='*80}")
        print("Signal Specificity Analysis:")
        
        if 'N-terminal' in region_accuracies:
            n_term_acc = region_accuracies['N-terminal']
            shuffled_acc = region_accuracies.get('Shuffled N-terminal', 0)
            c_term_acc = region_accuracies.get('C-terminal', 0)
            middle_acc = region_accuracies.get('Middle Region', 0)
            random_acc = region_accuracies.get('Random Fragment', 0)
            full_acc = region_accuracies.get('Full-length', 0)
            
            # Use multi-trial mean for random if available
            if 'Random Fragment (Multi-trial)' in region_accuracies:
                random_acc = region_accuracies['Random Fragment (Multi-trial)']
                print(f"Note: Using multi-trial mean for random fragment: {random_acc:.4f}")
            
            # Check if N-terminal significantly outperforms control regions
            n_vs_shuffled = n_term_acc - shuffled_acc
            n_vs_c = n_term_acc - c_term_acc
            n_vs_middle = n_term_acc - middle_acc
            n_vs_random = n_term_acc - random_acc
            
            print(f"1. N-terminal vs Shuffled N-terminal: {n_vs_shuffled:+.4f} (key comparison)")
            print(f"2. N-terminal vs C-terminal: {n_vs_c:+.4f}")
            print(f"3. N-terminal vs Middle Region: {n_vs_middle:+.4f}")
            print(f"4. N-terminal vs Random Fragment: {n_vs_random:+.4f}")
            print(f"5. N-terminal vs Full-length: {n_term_acc - full_acc:+.4f}")
            
            # Determine N-terminal signal specificity
            if (n_vs_shuffled > 0 and n_vs_c > 0 and n_vs_middle > 0 and n_vs_random > 0):
                print("\n✓ N-terminal signal shows strong specificity:")
                print("   - N-terminal accuracy higher than all control regions (including shuffled)")
                print("   - Indicates N-terminal contains specific sequence-order dependent classification information")
                print("   - Information is not just due to amino acid composition")
                
                if abs(n_term_acc - full_acc) < 0.05:
                    print(f"   - N-terminal performance comparable to full-length (difference: {abs(n_term_acc - full_acc):.4f})")
                    print(f"   - Indicates 100 N-terminal amino acids contain most classification information")
            elif n_vs_shuffled > 0:
                print("\n⚠️ N-terminal shows some sequence-order specificity:")
                print("   - N-terminal accuracy higher than shuffled N-terminal")
                print("   - Indicates sequence order (not just composition) contributes to classification")
            else:
                print("\n✗ N-terminal signal specificity insufficient:")
                print("   - N-terminal accuracy not higher than shuffled control")
                print("   - Classification may be based primarily on amino acid composition, not sequence order")
                print("   - N-terminal signal may not contain specific sequence-order information")
        
        # Performance ranking
        print(f"\n{'='*80}")
        print("Performance Ranking by Region (highest to lowest):")
        sorted_regions = sorted(region_accuracies.items(), key=lambda x: x[1], reverse=True)
        for i, (region, acc) in enumerate(sorted_regions, 1):
            print(f"{i}. {region}: {acc:.4f}")
        
        # Biological implications
        print(f"\n{'='*80}")
        print("Biological Implications for P450 Enzyme Functional Classification:")
        
        if 'N-terminal' in region_accuracies and region_accuracies['N-terminal'] == max(region_accuracies.values()):
            print("1. N-terminal amino acid sequences play key role in P450 enzyme functional classification")
            print("2. N-terminal contains sequence-order specific information, not just amino acid composition")
            print("3. The specific order of amino acids in N-terminal is important for functional classification")
            print("4. Targeted engineering of N-terminal may more effectively alter enzyme functional properties")
        elif 'Shuffled N-terminal' in region_accuracies and region_accuracies['Shuffled N-terminal'] == max(region_accuracies.values()):
            print("1. Amino acid composition (not sequence order) is sufficient for classification")
            print("2. N-terminal sequence order may not be critical for functional classification")
            print("3. Functional classification may depend more on overall amino acid composition")
        elif 'Full-length' in region_accuracies and region_accuracies['Full-length'] == max(region_accuracies.values()):
            print("1. Full-length sequences contain the most complete classification information")
            print("2. Functional classification information may be distributed across multiple domains")
            print("3. More comprehensive sequence analysis needed to understand functional classification mechanisms")
        else:
            print("1. Functional classification information may be distributed in specific sequence regions")
            print("2. Different regions may contain complementary classification information")
            print("3. Further analysis needed to understand structure-function relationships in specific regions")
    
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*80)
    print("Multi-region comparative analysis completed!")
    print(f"Detailed results available in directory: {OUTPUT_DIR}")
    print("="*80)

if __name__ == "__main__":
    main()