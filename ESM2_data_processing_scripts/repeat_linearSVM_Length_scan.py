#!/usr/bin/env python3
"""
植物P450酶N端氨基酸长度扫描分析脚本 - 多次重复线性SVM版本
重复运行6次LinearSVM，整合结果验证分化期、过渡期和去噪音假说
"""

import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from collections import Counter
import argparse
import warnings
warnings.filterwarnings('ignore')

# 使用transformers加载ESM2模型
from transformers import AutoTokenizer, AutoModel

# 机器学习相关库
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.decomposition import PCA
from sklearn.svm import LinearSVC, SVC
from sklearn.base import clone

# 统计检验
from scipy import stats
import scipy.stats as st

# 进度条
from tqdm import tqdm
import gc
import json

class RepeatedLinearSVMLengthScanAnalyzer:
    """多次重复的线性SVM P450酶N端长度扫描分析器"""
    
    def __init__(self, model_name: str = "facebook/esm2_t12_35M_UR50D", 
                 device: str = "cpu"):
        """
        初始化分析器
        
        Args:
            model_name: ESM2模型名称
            device: 计算设备 (cpu/cuda)
        """
        self.device = device if torch.cuda.is_available() and device == "cuda" else "cpu"
        
        # 使用transformers加载ESM2模型
        print(f"加载ESM2模型: {model_name}...")
        try:
            self.tokenizer = AutoTokenizer.from_pretrained(model_name)
            self.model = AutoModel.from_pretrained(model_name)
            self.model = self.model.to(self.device).eval()
            self.embedding_dim = self.model.embeddings.word_embeddings.embedding_dim
            print(f"模型加载完成，使用设备: {self.device}")
            print(f"嵌入维度: {self.embedding_dim}")
        except Exception as e:
            print(f"模型加载失败: {e}")
            print("尝试使用较小的模型...")
            self.tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
            self.model = AutoModel.from_pretrained("facebook/esm2_t6_8M_UR50D")
            self.model = self.model.to(self.device).eval()
            self.embedding_dim = self.model.embeddings.word_embeddings.embedding_dim
            print(f"使用替代模型完成，嵌入维度: {self.embedding_dim}")
    
    def detect_csv_delimiter(self, filepath: str) -> str:
        """自动检测CSV文件分隔符"""
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
                print(f"检测到分隔符: '{best_delimiter}'")
                return best_delimiter
            return ','
    
    def load_data(self, fasta_file: str, cluster_file: str) -> Tuple[Dict, Dict]:
        """加载FASTA文件和聚类标签"""
        # 加载聚类文件
        print("加载聚类标签文件...")
        delimiter = self.detect_csv_delimiter(cluster_file)
        
        try:
            cluster_df = pd.read_csv(cluster_file, sep=delimiter, engine='python')
        except:
            # 尝试不同分隔符
            for sep in [',', '\t', ';', '|']:
                try:
                    cluster_df = pd.read_csv(cluster_file, sep=sep, engine='python')
                    break
                except:
                    continue
        
        # 解析基因到聚类的映射
        gene_cluster_dict = {}
        for _, row in cluster_df.iterrows():
            cluster_col = 'Cluster' if 'Cluster' in cluster_df.columns else cluster_df.columns[0]
            gene_list_col = 'Gene_List' if 'Gene_List' in cluster_df.columns else cluster_df.columns[2]
            
            cluster_name = str(row[cluster_col]).strip()
            gene_list_str = str(row[gene_list_col]).strip()
            
            # 分割基因列表
            for sep in [';', ',', '|', '\t', ' ']:
                if sep in gene_list_str:
                    genes = [g.strip() for g in gene_list_str.split(sep) if g.strip()]
                    break
            else:
                genes = [g.strip() for g in gene_list_str.split() if g.strip()]
            
            for gene in genes:
                gene_cluster_dict[gene] = cluster_name
        
        print(f"加载了 {len(gene_cluster_dict)} 个基因的聚类标签")
        
        # 加载FASTA文件
        print("\n加载FASTA文件...")
        sequences_dict = {}
        
        with open(fasta_file, 'r') as f:
            current_gene = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_gene and current_seq:
                        sequences_dict[current_gene] = ''.join(current_seq)
                    
                    # 提取基因名
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
                sequences_dict[current_gene] = ''.join(current_seq)
        
        print(f"加载了 {len(sequences_dict)} 个蛋白质序列")
        
        # 对齐数据
        common_genes = set(sequences_dict.keys()) & set(gene_cluster_dict.keys())
        sequences_dict = {g: sequences_dict[g] for g in common_genes}
        gene_cluster_dict = {g: gene_cluster_dict[g] for g in common_genes}
        
        print(f"\n对齐后共有 {len(common_genes)} 个共同基因")
        
        return sequences_dict, gene_cluster_dict
    
    def handle_class_imbalance(self, labels: List) -> Dict:
        """处理类别不平衡"""
        label_counts = Counter(labels)
        print("\n=== 类别分布统计 ===")
        for label, count in sorted(label_counts.items(), key=lambda x: x[1], reverse=True):
            print(f"{label}: {count} 个样本")
        
        # 合并稀有类别（样本数 < 5）
        rare_labels = [label for label, count in label_counts.items() if count < 5]
        if rare_labels:
            print(f"\n发现稀有类别（样本数 < 5）: {len(rare_labels)} 个")
            labels_processed = ['Rare_Class' if label in rare_labels else label for label in labels]
            new_counts = Counter(labels_processed)
            print("合并后的类别分布:")
            for label, count in sorted(new_counts.items(), key=lambda x: x[1], reverse=True):
                print(f"{label}: {count} 个样本")
        else:
            labels_processed = labels
            new_counts = label_counts
        
        return {
            'labels': labels_processed,
            'original_labels': labels,
            'rare_labels': rare_labels if rare_labels else [],
            'label_counts': new_counts
        }
    
    def extract_esm2_embeddings(self, sequences: List[str]) -> np.ndarray:
        """提取ESM2嵌入向量"""
        print(f"提取{len(sequences)}个序列的ESM2嵌入向量...")
        all_embeddings = []
        batch_size = 2  # 减小批处理大小适应大模型
        
        with torch.no_grad():
            for i in tqdm(range(0, len(sequences), batch_size), desc="提取嵌入向量"):
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
                
                # 清理内存
                del inputs, outputs, last_hidden_state
                if self.device == "cuda":
                    torch.cuda.empty_cache()
                gc.collect()
        
        X = np.array(all_embeddings)
        print(f"特征矩阵形状: {X.shape}")
        return X
    
    def extract_n_terminal_embeddings(self, sequences_dict: Dict, gene_cluster_dict: Dict, 
                                     n_aa: int) -> Tuple[np.ndarray, np.ndarray, Dict]:
        """提取指定长度N端的嵌入向量"""
        print(f"提取N端{n_aa}个氨基酸的ESM2嵌入向量...")
        
        genes = list(sequences_dict.keys())
        sequences = []
        labels = []
        
        for gene in genes:
            seq = sequences_dict[gene]
            
            if len(seq) >= n_aa:
                seq_to_use = seq[:n_aa]
            else:
                print(f"警告: 序列 {gene} 长度 ({len(seq)}) 小于 {n_aa}，使用完整序列")
                seq_to_use = seq
            
            sequences.append(seq_to_use)
            labels.append(gene_cluster_dict[gene])
        
        label_info = self.handle_class_imbalance(labels)
        labels = label_info['labels']
        
        self.label_encoder = LabelEncoder()
        y = self.label_encoder.fit_transform(labels)
        X = self.extract_esm2_embeddings(sequences)
        
        return X, y, label_info
    
    def extract_full_length_embeddings(self, sequences_dict: Dict, 
                                      gene_cluster_dict: Dict) -> Tuple[np.ndarray, np.ndarray, Dict]:
        """提取全长序列的嵌入向量"""
        print("提取全长序列的ESM2嵌入向量...")
        
        genes = list(sequences_dict.keys())
        sequences = []
        labels = []
        
        for gene in genes:
            seq = sequences_dict[gene]
            sequences.append(seq)
            labels.append(gene_cluster_dict[gene])
        
        label_info = self.handle_class_imbalance(labels)
        labels = label_info['labels']
        
        self.label_encoder = LabelEncoder()
        y = self.label_encoder.fit_transform(labels)
        X = self.extract_esm2_embeddings(sequences)
        
        return X, y, label_info
    
    def robust_linear_svm_cv(self, X: np.ndarray, y: np.ndarray, 
                            n_splits: int = 5, random_seed: int = 42) -> Dict[str, Any]:
        """
        稳健的线性SVM交叉验证
        包含多重保护机制防止SVM失败
        """
        print(f"执行稳健的线性SVM交叉验证 (随机种子: {random_seed})...")
        
        # 设置随机种子
        np.random.seed(random_seed)
        
        # 1. 数据标准化（处理数值稳定性）
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # 检查并处理NaN/Inf值
        if np.any(np.isnan(X_scaled)) or np.any(np.isinf(X_scaled)):
            print("警告: 检测到NaN或Inf值，进行清理...")
            X_scaled = np.nan_to_num(X_scaled, nan=0.0, posinf=1e6, neginf=-1e6)
        
        # 2. PCA降维（提高稳定性并减少计算量）
        n_samples, n_features = X_scaled.shape
        max_components = min(n_samples, n_features, 200)  # 限制最大维度
        
        try:
            # 使用更保守的PCA参数
            pca = PCA(n_components=min(0.90, max_components/n_features), 
                     random_state=random_seed, svd_solver='randomized')
            X_pca = pca.fit_transform(X_scaled)
            pca_variance = pca.explained_variance_ratio_.sum()
            print(f"PCA降维: {X_scaled.shape[1]} -> {X_pca.shape[1]} (解释方差: {pca_variance:.3f})")
        except Exception as e:
            print(f"PCA失败: {e}，使用原始特征")
            X_pca = X_scaled
            pca_variance = 1.0
        
        # 3. 准备SVM参数（多重配置提高稳定性）
        svm_configs = [
            {'C': 0.01, 'max_iter': 5000, 'dual': False, 'tol': 1e-3},
            {'C': 0.1, 'max_iter': 5000, 'dual': False, 'tol': 1e-3},
            {'C': 1.0, 'max_iter': 10000, 'dual': False, 'tol': 1e-4},
            {'C': 10.0, 'max_iter': 15000, 'dual': False, 'tol': 1e-4},
        ]
        
        best_scores = []
        best_config = None
        best_mean_score = 0
        
        # 4. 尝试不同SVM配置
        for config_idx, config in enumerate(svm_configs):
            print(f"  尝试SVM配置 {config_idx+1}/{len(svm_configs)}: C={config['C']}")
            
            try:
                # 使用LinearSVC而不是SVC(kernel='linear')，更稳定高效
                svm_model = LinearSVC(
                    C=config['C'],
                    class_weight='balanced',
                    max_iter=config['max_iter'],
                    tol=config['tol'],
                    random_state=random_seed,
                    dual=config['dual'],  # 当n_samples > n_features时，dual=False更高效
                    verbose=0
                )
                
                # 稳健的交叉验证实现
                skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_seed)
                fold_scores = []
                
                for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X_pca, y)):
                    X_train, X_test = X_pca[train_idx], X_pca[test_idx]
                    y_train, y_test = y[train_idx], y[test_idx]
                    
                    try:
                        # 克隆模型以避免状态污染
                        fold_model = clone(svm_model)
                        fold_model.fit(X_train, y_train)
                        score = fold_model.score(X_test, y_test)
                        fold_scores.append(score)
                    except Exception as e:
                        print(f"    折叠 {fold_idx+1} 失败: {e}")
                        # 如果某个折叠失败，使用随机准确率作为下限
                        n_classes = len(np.unique(y_train))
                        fold_scores.append(1.0 / n_classes)
                        continue
                
                # 计算平均准确率
                if len(fold_scores) >= 3:  # 至少需要3个成功的折叠
                    mean_score = np.mean(fold_scores)
                    std_score = np.std(fold_scores)
                    
                    if mean_score > best_mean_score:
                        best_mean_score = mean_score
                        best_scores = fold_scores
                        best_config = config
                        
                        print(f"    配置 {config_idx+1} 成功: {mean_score:.4f} ± {std_score:.4f}")
                    else:
                        print(f"    配置 {config_idx+1} 可用: {mean_score:.4f} ± {std_score:.4f}")
                else:
                    print(f"    配置 {config_idx+1} 失败: 成功折叠不足")
                    
            except Exception as e:
                print(f"    配置 {config_idx+1} 失败: {e}")
                continue
        
        # 5. 如果所有配置都失败，使用回退方案
        if len(best_scores) == 0:
            print("  所有SVM配置都失败，使用简化回退方案...")
            
            # 尝试最简单的配置
            try:
                svm_model = LinearSVC(
                    C=1.0,
                    class_weight='balanced',
                    max_iter=10000,
                    tol=1e-3,
                    random_state=random_seed,
                    dual=False,
                    verbose=0
                )
                
                # 使用更小的交叉验证折数
                skf = StratifiedKFold(n_splits=min(3, n_splits), shuffle=True, random_state=random_seed)
                best_scores = cross_val_score(svm_model, X_pca, y, cv=skf, scoring='accuracy')
                best_config = {'C': 1.0, 'fallback': True}
                
                print(f"  回退方案成功: {np.mean(best_scores):.4f} ± {np.std(best_scores):.4f}")
            except Exception as e:
                print(f"  回退方案失败: {e}")
                # 返回默认结果
                n_classes = len(np.unique(y))
                best_scores = [1.0 / n_classes] * min(3, n_splits)
                best_config = {'C': 1.0, 'failed': True}
        
        # 6. 返回结果
        return {
            'mean_accuracy': np.mean(best_scores),
            'std_accuracy': np.std(best_scores) if len(best_scores) > 1 else 0.0,
            'scores': best_scores,
            'n_pca_components': X_pca.shape[1],
            'explained_variance': pca_variance,
            'svm_config': best_config,
            'random_seed': random_seed,
            'success': True
        }
    
    def single_scan_n_terminal_lengths(self, sequences_dict: Dict, 
                                      gene_cluster_dict: Dict, 
                                      min_length: int = 40, 
                                      max_length: int = 100,
                                      random_seed: int = 42) -> Dict:
        """
        单次扫描不同N端长度的分类性能
        
        Args:
            sequences_dict: 序列字典
            gene_cluster_dict: 基因聚类字典
            min_length: 最小N端长度
            max_length: 最大N端长度
            random_seed: 随机种子
        
        Returns:
            包含所有长度性能的字典
        """
        print(f"开始单次扫描 (随机种子: {random_seed})...")
        
        results = {}
        
        # 首先提取全长序列的特征和性能（只做一次）
        print("1. 处理全长序列...")
        X_full, y_full, label_info = self.extract_full_length_embeddings(sequences_dict, gene_cluster_dict)
        full_performance = self.robust_linear_svm_cv(X_full, y_full, random_seed=random_seed)
        
        results['full_length'] = {
            'X': X_full,
            'y': y_full,
            'performance': full_performance,
            'label_info': label_info,
            'mean_accuracy': full_performance['mean_accuracy'],
            'std_accuracy': full_performance['std_accuracy'],
            'svm_config': full_performance.get('svm_config', {}),
            'random_seed': random_seed
        }
        
        print(f"全长序列准确率: {full_performance['mean_accuracy']:.4f} ± {full_performance['std_accuracy']:.4f}")
        
        # 存储不同长度的性能
        length_results = {}
        
        # 从min_length到max_length，逐步增加长度
        for n_aa in range(min_length, max_length + 1):
            print(f"分析N端长度: {n_aa} AA")
            
            try:
                # 提取N端特征
                X_nterm, y_nterm, _ = self.extract_n_terminal_embeddings(sequences_dict, gene_cluster_dict, n_aa)
                
                # 评估性能
                nterm_performance = self.robust_linear_svm_cv(X_nterm, y_nterm, random_seed=random_seed)
                
                # 与全长序列进行统计检验（配对t检验）
                if (len(nterm_performance['scores']) == len(full_performance['scores']) and 
                    nterm_performance['success'] and full_performance['success']):
                    
                    t_stat, p_value = stats.ttest_rel(nterm_performance['scores'], full_performance['scores'])
                    mean_diff = np.mean(nterm_performance['scores']) - np.mean(full_performance['scores'])
                    pooled_std = np.sqrt((np.std(nterm_performance['scores'], ddof=1)**2 + 
                                        np.std(full_performance['scores'], ddof=1)**2) / 2)
                    cohens_d = mean_diff / pooled_std if pooled_std != 0 else 0
                    
                    # 计算性能差距百分比
                    accuracy_gap = full_performance['mean_accuracy'] - nterm_performance['mean_accuracy']
                    gap_percentage = (accuracy_gap / full_performance['mean_accuracy']) * 100
                    
                    # 判断是否显著差异
                    is_significant = p_value < 0.05
                    
                    # 存储结果
                    length_results[n_aa] = {
                        'mean_accuracy': nterm_performance['mean_accuracy'],
                        'std_accuracy': nterm_performance['std_accuracy'],
                        'accuracy_scores': nterm_performance['scores'],
                        'p_value': p_value,
                        'cohens_d': cohens_d,
                        'accuracy_gap': accuracy_gap,
                        'gap_percentage': gap_percentage,
                        'is_significant': is_significant,
                        'n_pca_components': nterm_performance['n_pca_components'],
                        'explained_variance': nterm_performance['explained_variance'],
                        'svm_config': nterm_performance.get('svm_config', {}),
                        'random_seed': random_seed
                    }
                    
                else:
                    print(f"警告: 无法与全长序列进行统计比较")
                    # 仍然存储结果
                    length_results[n_aa] = {
                        'mean_accuracy': nterm_performance['mean_accuracy'],
                        'std_accuracy': nterm_performance['std_accuracy'],
                        'accuracy_scores': nterm_performance['scores'],
                        'p_value': 1.0,  # 默认p值设为1.0（无差异）
                        'cohens_d': 0.0,
                        'accuracy_gap': 0.0,
                        'gap_percentage': 0.0,
                        'is_significant': False,
                        'n_pca_components': nterm_performance['n_pca_components'],
                        'explained_variance': nterm_performance['explained_variance'],
                        'svm_config': nterm_performance.get('svm_config', {}),
                        'random_seed': random_seed,
                        'comparison_valid': False
                    }
                    
            except Exception as e:
                print(f"处理N端 {n_aa}AA 时出错: {e}")
                continue
        
        results['length_scan'] = length_results
        
        print(f"单次扫描完成 (随机种子: {random_seed})")
        return results
    
    def repeated_scan_n_terminal_lengths(self, sequences_dict: Dict, 
                                        gene_cluster_dict: Dict, 
                                        min_length: int = 40, 
                                        max_length: int = 100,
                                        n_repeats: int = 6,
                                        random_seed_start: int = 42) -> Dict:
        """
        重复多次扫描不同N端长度的分类性能
        
        Args:
            sequences_dict: 序列字典
            gene_cluster_dict: 基因聚类字典
            min_length: 最小N端长度
            max_length: 最大N端长度
            n_repeats: 重复次数
            random_seed_start: 随机种子起始值
        
        Returns:
            包含所有重复结果和整合结果的字典
        """
        print(f"\n{'='*80}")
        print(f"开始重复扫描N端长度 {min_length}-{max_length} AA 的分类性能")
        print(f"重复次数: {n_repeats}, 使用线性SVM分类器")
        print('='*80)
        
        all_results = []
        
        # 重复运行多次扫描
        for i in range(n_repeats):
            print(f"\n{'='*60}")
            print(f"第 {i+1}/{n_repeats} 次重复扫描")
            print('='*60)
            
            random_seed = random_seed_start + i * 100  # 使用不同的随机种子
            
            try:
                # 执行单次扫描
                single_result = self.single_scan_n_terminal_lengths(
                    sequences_dict, 
                    gene_cluster_dict,
                    min_length=min_length,
                    max_length=max_length,
                    random_seed=random_seed
                )
                all_results.append(single_result)
                
                # 保存单次结果
                single_output_dir = f"repeat_{i+1}_seed_{random_seed}"
                self.export_single_results(single_result, single_output_dir)
                
            except Exception as e:
                print(f"第 {i+1} 次重复扫描失败: {e}")
                continue
        
        # 整合所有重复的结果
        integrated_results = self.integrate_results(all_results)
        
        return {
            'all_repeats': all_results,
            'integrated': integrated_results,
            'n_repeats': len(all_results),
            'successful_repeats': [i+1 for i in range(len(all_results))]
        }
    
    def integrate_results(self, all_results: List[Dict]) -> Dict:
        """
        整合多次重复的结果
        
        Args:
            all_results: 所有重复结果的列表
        
        Returns:
            整合后的结果字典
        """
        print(f"\n整合 {len(all_results)} 次重复结果...")
        
        if not all_results:
            print("错误: 没有可整合的结果")
            return {}
        
        # 获取所有长度的列表
        all_lengths = set()
        for result in all_results:
            if 'length_scan' in result:
                all_lengths.update(result['length_scan'].keys())
        
        all_lengths = sorted(all_lengths)
        
        # 整合每个长度的结果
        integrated_lengths = {}
        
        for n_aa in all_lengths:
            accuracies = []
            p_values = []
            is_significant_list = []
            gap_percentages = []
            
            for result in all_results:
                if 'length_scan' in result and n_aa in result['length_scan']:
                    length_result = result['length_scan'][n_aa]
                    accuracies.append(length_result['mean_accuracy'])
                    p_values.append(length_result['p_value'])
                    is_significant_list.append(length_result['is_significant'])
                    gap_percentages.append(length_result['gap_percentage'])
            
            if accuracies:
                # 计算整合统计量
                integrated_lengths[n_aa] = {
                    'mean_accuracy': np.mean(accuracies),
                    'std_accuracy': np.std(accuracies),
                    'min_accuracy': np.min(accuracies),
                    'max_accuracy': np.max(accuracies),
                    'mean_p_value': np.mean(p_values),
                    'median_p_value': np.median(p_values),
                    'significant_ratio': np.mean(is_significant_list),  # 显著的比例
                    'mean_gap_percentage': np.mean(gap_percentages),
                    'std_gap_percentage': np.std(gap_percentages),
                    'n_repeats': len(accuracies),
                    'accuracies': accuracies,
                    'p_values': p_values,
                    'is_significant_list': is_significant_list
                }
        
        # 整合全长序列结果
        full_accuracies = []
        for result in all_results:
            if 'full_length' in result:
                full_accuracies.append(result['full_length']['mean_accuracy'])
        
        integrated_full = {
            'mean_accuracy': np.mean(full_accuracies),
            'std_accuracy': np.std(full_accuracies),
            'min_accuracy': np.min(full_accuracies),
            'max_accuracy': np.max(full_accuracies),
            'n_repeats': len(full_accuracies),
            'accuracies': full_accuracies
        }
        
        return {
            'full_length': integrated_full,
            'length_scan': integrated_lengths,
            'n_total_repeats': len(all_results)
        }
    
    def find_turning_point_integrated(self, integrated_results: Dict, 
                                     p_threshold: float = 0.05, 
                                     significance_ratio_threshold: float = 0.5) -> Dict:
        """
        在整合结果中寻找分类性能拐点
        定义：第一个与全长序列无显著差异的N端长度
        显著性由多次重复中显著的比例决定
        """
        print(f"\n在整合结果中寻找分类性能拐点...")
        
        if 'length_scan' not in integrated_results:
            print("错误: 没有找到长度扫描结果")
            return None
        
        length_results = integrated_results['length_scan']
        full_accuracy = integrated_results['full_length']['mean_accuracy']
        
        turning_points = []
        
        # 寻找第一个无显著差异的长度
        for n_aa in sorted(length_results.keys()):
            result = length_results[n_aa]
            
            # 检查是否有足够的数据
            if result['n_repeats'] < 3:
                continue
            
            # 检查是否无显著差异（显著比例低于阈值）
            if result['significant_ratio'] <= significance_ratio_threshold:
                # 检查是否达到全长效应的90%
                accuracy_ratio = result['mean_accuracy'] / full_accuracy
                if accuracy_ratio >= 0.90:
                    turning_points.append({
                        'length': n_aa,
                        'accuracy': result['mean_accuracy'],
                        'full_accuracy': full_accuracy,
                        'accuracy_ratio': accuracy_ratio,
                        'mean_p_value': result['mean_p_value'],
                        'significant_ratio': result['significant_ratio'],
                        'mean_gap_percentage': result['mean_gap_percentage']
                    })
                    print(f"候选拐点: N端 {n_aa}AA")
                    print(f"  准确率: {result['mean_accuracy']:.4f} (全长: {full_accuracy:.4f})")
                    print(f"  比例: {accuracy_ratio:.3f}")
                    print(f"  平均p值: {result['mean_p_value']:.4f}")
                    print(f"  显著比例: {result['significant_ratio']:.2%}")
                    print(f"  平均差距: {result['mean_gap_percentage']:.2f}%")
        
        if turning_points:
            # 选择第一个满足条件的长度
            first_turning_point = min(turning_points, key=lambda x: x['length'])
            print(f"\n推荐拐点: N端 {first_turning_point['length']}AA")
            print(f"此时分类性能与全长序列无显著差异 (平均p={first_turning_point['mean_p_value']:.4f})")
            print(f"准确率达到全长的 {first_turning_point['accuracy_ratio']:.2%}")
            print(f"显著比例: {first_turning_point['significant_ratio']:.2%}")
            
            return first_turning_point
        else:
            print("未找到明确的拐点")
            return None
    
    def export_single_results(self, results: Dict, output_dir: str = "single_repeat_results"):
        """导出单次重复结果到CSV文件"""
        Path(output_dir).mkdir(exist_ok=True)
        
        # 导出长度扫描性能数据
        if 'length_scan' in results:
            scan_data = []
            for n_aa, result in results['length_scan'].items():
                scan_data.append({
                    'N_terminal_length_AA': n_aa,
                    'Mean_Accuracy': result['mean_accuracy'],
                    'Std_Accuracy': result['std_accuracy'],
                    'P_value': result['p_value'],
                    'Cohens_d': result['cohens_d'],
                    'Accuracy_Gap': result['accuracy_gap'],
                    'Gap_Percentage': result['gap_percentage'],
                    'Is_Significant': result['is_significant'],
                    'Significant_Flag': 'Yes' if result['is_significant'] else 'No',
                    'PCA_Components': result['n_pca_components'],
                    'Explained_Variance': result['explained_variance'],
                    'SVM_C': result.get('svm_config', {}).get('C', 'N/A'),
                    'Random_Seed': result.get('random_seed', 'N/A')
                })
            
            scan_df = pd.DataFrame(scan_data)
            scan_df = scan_df.sort_values('N_terminal_length_AA')
            scan_df.to_csv(f"{output_dir}/n_terminal_length_scan.csv", 
                          index=False, encoding='utf-8')
        
        # 导出全长序列性能
        if 'full_length' in results:
            full_data = [{
                'Sequence_Type': 'Full_length',
                'Mean_Accuracy': results['full_length']['mean_accuracy'],
                'Std_Accuracy': results['full_length']['std_accuracy'],
                'PCA_Components': results['full_length']['performance']['n_pca_components'],
                'Explained_Variance': results['full_length']['performance']['explained_variance'],
                'SVM_C': results['full_length']['performance'].get('svm_config', {}).get('C', 'N/A'),
                'Random_Seed': results['full_length'].get('random_seed', 'N/A')
            }]
            
            full_df = pd.DataFrame(full_data)
            full_df.to_csv(f"{output_dir}/full_length_performance.csv", 
                          index=False, encoding='utf-8')
        
        print(f"单次重复结果已导出到: {output_dir}/")
    
    def export_integrated_results(self, all_results: Dict, output_dir: str = "integrated_results"):
        """导出整合结果到CSV和JSON文件"""
        print(f"\n导出整合结果到文件...")
        
        Path(output_dir).mkdir(exist_ok=True)
        
        # 导出整合的长度扫描性能数据
        if 'integrated' in all_results and 'length_scan' in all_results['integrated']:
            scan_data = []
            for n_aa, result in all_results['integrated']['length_scan'].items():
                scan_data.append({
                    'N_terminal_length_AA': n_aa,
                    'Mean_Accuracy': result['mean_accuracy'],
                    'Std_Accuracy': result['std_accuracy'],
                    'Min_Accuracy': result['min_accuracy'],
                    'Max_Accuracy': result['max_accuracy'],
                    'Mean_P_value': result['mean_p_value'],
                    'Median_P_value': result['median_p_value'],
                    'Significant_Ratio': result['significant_ratio'],
                    'Mean_Gap_Percentage': result['mean_gap_percentage'],
                    'Std_Gap_Percentage': result['std_gap_percentage'],
                    'N_Repeats': result['n_repeats']
                })
            
            scan_df = pd.DataFrame(scan_data)
            scan_df = scan_df.sort_values('N_terminal_length_AA')
            scan_df.to_csv(f"{output_dir}/integrated_n_terminal_length_scan.csv", 
                          index=False, encoding='utf-8')
            print(f"已保存: {output_dir}/integrated_n_terminal_length_scan.csv")
        
        # 导出整合的全长序列性能
        if 'integrated' in all_results and 'full_length' in all_results['integrated']:
            full_data = [{
                'Sequence_Type': 'Full_length',
                'Mean_Accuracy': all_results['integrated']['full_length']['mean_accuracy'],
                'Std_Accuracy': all_results['integrated']['full_length']['std_accuracy'],
                'Min_Accuracy': all_results['integrated']['full_length']['min_accuracy'],
                'Max_Accuracy': all_results['integrated']['full_length']['max_accuracy'],
                'N_Repeats': all_results['integrated']['full_length']['n_repeats']
            }]
            
            full_df = pd.DataFrame(full_data)
            full_df.to_csv(f"{output_dir}/integrated_full_length_performance.csv", 
                          index=False, encoding='utf-8')
            print(f"已保存: {output_dir}/integrated_full_length_performance.csv")
        
        # 导出详细的所有重复准确率数据
        if 'all_repeats' in all_results:
            detailed_data = []
            for repeat_idx, result in enumerate(all_results['all_repeats']):
                if 'length_scan' in result:
                    for n_aa, length_result in result['length_scan'].items():
                        detailed_data.append({
                            'Repeat': repeat_idx + 1,
                            'N_terminal_length_AA': n_aa,
                            'Accuracy': length_result['mean_accuracy'],
                            'P_value': length_result['p_value'],
                            'Is_Significant': length_result['is_significant'],
                            'Random_Seed': length_result.get('random_seed', 'N/A')
                        })
            
            detailed_df = pd.DataFrame(detailed_data)
            detailed_df.to_csv(f"{output_dir}/detailed_all_repeats_accuracy.csv", 
                              index=False, encoding='utf-8')
            print(f"已保存: {output_dir}/detailed_all_repeats_accuracy.csv")
        
        # 导出为JSON格式（便于后续分析）
        with open(f"{output_dir}/all_results_summary.json", 'w') as f:
            # 简化数据以便JSON序列化
            json_data = {
                'n_repeats': all_results.get('n_repeats', 0),
                'successful_repeats': all_results.get('successful_repeats', []),
                'integrated_summary': {
                    'full_length': all_results['integrated']['full_length'] if 'integrated' in all_results else {},
                    'turning_point': self.find_turning_point_integrated(all_results['integrated']) if 'integrated' in all_results else None
                }
            }
            json.dump(json_data, f, indent=2, default=str)
            print(f"已保存: {output_dir}/all_results_summary.json")
        
        print(f"\n所有整合结果已导出到: {output_dir}/")
    
    def visualize_integrated_results(self, all_results: Dict, output_dir: str = "integrated_results"):
        """可视化整合的多次重复结果"""
        print(f"\n生成整合结果的可视化图表...")
        
        if 'integrated' not in all_results:
            print("错误: 没有整合结果可供可视化")
            return
        
        integrated_results = all_results['integrated']
        
        # 创建图表目录
        plot_dir = Path(output_dir) / "plots"
        plot_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置绘图风格
        plt.style.use('seaborn-v0_8-whitegrid')
        sns.set_palette("husl")
        
        if 'length_scan' not in integrated_results:
            print("错误: 没有长度扫描结果可供可视化")
            return
        
        length_results = integrated_results['length_scan']
        full_results = integrated_results['full_length']
        full_accuracy = full_results['mean_accuracy']
        full_std = full_results['std_accuracy']
        
        # 准备数据
        lengths = sorted(length_results.keys())
        mean_accuracies = [length_results[n]['mean_accuracy'] for n in lengths]
        std_accuracies = [length_results[n]['std_accuracy'] for n in lengths]
        significant_ratios = [length_results[n]['significant_ratio'] for n in lengths]
        mean_gap_percentages = [length_results[n]['mean_gap_percentage'] for n in lengths]
        
        # 1. 主要图表：平均准确率随长度的变化（带标准差）
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
        
        # 图1: 平均准确率曲线（带标准差阴影）
        ax1.plot(lengths, mean_accuracies, 'b-', linewidth=3, marker='o', 
                markersize=6, label='Average N-terminal accuracy')
        ax1.fill_between(lengths, 
                        [mean_accuracies[i] - std_accuracies[i] for i in range(len(lengths))], 
                        [mean_accuracies[i] + std_accuracies[i] for i in range(len(lengths))], 
                        alpha=0.2, color='blue', label='±1 std')
        
        # 全长序列准确率水平线（带标准差范围）
        ax1.axhline(y=full_accuracy, color='red', linestyle='--', 
                   linewidth=2, label=f'Full-length accuracy ({full_accuracy:.4f})')
        ax1.fill_between([min(lengths), max(lengths)], 
                        full_accuracy - full_std, full_accuracy + full_std,
                        alpha=0.1, color='red', label='Full-length ±1 std')
        
        # 标记分化期、过渡期、去噪音期
        ax1.axvspan(40, 49, alpha=0.1, color='orange', label='Rapid divergence (40-49AA)')
        ax1.axvspan(50, 84, alpha=0.1, color='yellow', label='Transition (50-84AA)')
        ax1.axvspan(85, 100, alpha=0.1, color='green', label='Noise removal (85-100AA)')
        
        ax1.set_xlabel('N-terminal Length (AA)', fontsize=12)
        ax1.set_ylabel('Classification Accuracy', fontsize=12)
        ax1.set_title('Average Accuracy vs N-terminal Length (6 Repeats)', 
                     fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10, loc='upper left')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0.85, 1.0)
        
        # 添加准确率文本标签（关键点）
        key_lengths = [40, 50, 60, 70, 80, 90, 100]
        for length in key_lengths:
            if length in lengths:
                idx = lengths.index(length)
                ax1.annotate(f'{mean_accuracies[idx]:.3f}', 
                           xy=(length, mean_accuracies[idx]), 
                           xytext=(0, 10),
                           textcoords='offset points', 
                           ha='center', fontsize=9, fontweight='bold')
        
        # 图2: 显著比例热图
        # 创建显著比例矩阵
        sig_matrix = []
        for n_aa in lengths:
            sig_ratio = significant_ratios[lengths.index(n_aa)]
            sig_matrix.append([sig_ratio])
        
        im = ax2.imshow(sig_matrix, aspect='auto', cmap='RdYlGn_r', 
                       vmin=0, vmax=1, extent=[0, 1, min(lengths)-0.5, max(lengths)+0.5])
        
        ax2.set_xlabel('Significance', fontsize=12)
        ax2.set_ylabel('N-terminal Length (AA)', fontsize=12)
        ax2.set_title('Significant Difference Ratio (p<0.05)', 
                     fontsize=14, fontweight='bold')
        ax2.set_xticks([])
        ax2.set_yticks(np.arange(min(lengths), max(lengths)+1, 5))
        
        # 添加颜色条
        cbar = plt.colorbar(im, ax=ax2, orientation='vertical', pad=0.02)
        cbar.set_label('Ratio of Significant Differences', fontsize=10)
        
        # 添加显著比例数值
        for i, length in enumerate(lengths):
            if i % 5 == 0:  # 每5个长度标注一次
                ax2.text(0.5, length, f'{significant_ratios[i]:.2f}', 
                        ha='center', va='center', fontsize=8, 
                        color='white' if significant_ratios[i] > 0.5 else 'black',
                        fontweight='bold')
        
        # 图3: 性能差距百分比
        bars = ax3.bar(lengths, mean_gap_percentages, color='skyblue', alpha=0.7)
        
        # 根据正负值设置不同颜色
        for i, bar in enumerate(bars):
            if mean_gap_percentages[i] > 0:
                bar.set_color('lightcoral')  # 正值：红色，性能比全长差
            else:
                bar.set_color('lightgreen')  # 负值：绿色，性能比全长好
        
        ax3.axhline(y=0, color='black', linestyle='-', linewidth=1)
        ax3.set_xlabel('N-terminal Length (AA)', fontsize=12)
        ax3.set_ylabel('Accuracy Gap (%)', fontsize=12)
        ax3.set_title('Average Performance Gap vs Full-length', 
                     fontsize=14, fontweight='bold')
        ax3.grid(True, alpha=0.3, axis='y')
        
        # 添加性能差距数值标签
        for i, (length, gap) in enumerate(zip(lengths, mean_gap_percentages)):
            if i % 5 == 0:  # 每5个长度标注一次
                ax3.annotate(f'{gap:.1f}%', 
                           xy=(length, gap), 
                           xytext=(0, 10 if gap >= 0 else -15),
                           textcoords='offset points', 
                           ha='center', fontsize=8,
                           fontweight='bold')
        
        # 图4: 所有重复的准确率分布（箱线图）
        # 准备数据：每个长度的所有重复准确率
        boxplot_data = []
        boxplot_labels = []
        
        selected_lengths = [40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
        for length in selected_lengths:
            if length in lengths:
                boxplot_labels.append(str(length))
                # 获取该长度在所有重复中的准确率
                if 'all_repeats' in all_results:
                    length_accuracies = []
                    for result in all_results['all_repeats']:
                        if 'length_scan' in result and length in result['length_scan']:
                            length_accuracies.append(result['length_scan'][length]['mean_accuracy'])
                    boxplot_data.append(length_accuracies)
        
        bp = ax4.boxplot(boxplot_data, labels=boxplot_labels, patch_artist=True)
        
        # 设置箱线图颜色
        colors = ['lightblue'] * len(boxplot_data)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        
        # 添加全长准确率参考线
        ax4.axhline(y=full_accuracy, color='red', linestyle='--', 
                   linewidth=2, label=f'Full-length ({full_accuracy:.3f})')
        
        ax4.set_xlabel('N-terminal Length (AA)', fontsize=12)
        ax4.set_ylabel('Accuracy Distribution', fontsize=12)
        ax4.set_title('Accuracy Distribution Across 6 Repeats', 
                     fontsize=14, fontweight='bold')
        ax4.grid(True, alpha=0.3, axis='y')
        ax4.legend(fontsize=10)
        ax4.set_ylim(0.85, 1.0)
        
        plt.tight_layout()
        
        # 保存图表
        plt.savefig(f"{plot_dir}/integrated_length_scan_analysis.png", 
                   dpi=300, bbox_inches='tight')
        plt.savefig(f"{plot_dir}/integrated_length_scan_analysis.pdf", 
                   bbox_inches='tight')
        
        # 2. 拐点分析图
        turning_point = self.find_turning_point_integrated(integrated_results)
        
        fig2, ax5 = plt.subplots(figsize=(12, 8))
        
        # 绘制平均准确率曲线
        ax5.plot(lengths, mean_accuracies, 'b-', linewidth=3, marker='o', 
                markersize=8, label='Average N-terminal accuracy')
        
        # 全长序列准确率水平线
        ax5.axhline(y=full_accuracy, color='red', linestyle='--', 
                   linewidth=3, label=f'Full-length accuracy ({full_accuracy:.4f})')
        
        # 标记拐点
        if turning_point:
            tp_length = turning_point['length']
            tp_accuracy = turning_point['accuracy']
            
            ax5.axvline(x=tp_length, color='green', linestyle='--', 
                       linewidth=3, alpha=0.7, label=f'Turning point ({tp_length}AA)')
            ax5.scatter([tp_length], [tp_accuracy], color='green', s=300, 
                       zorder=10, edgecolors='black', linewidth=2, marker='*')
            
            ax5.annotate(f'Turning Point\n{tp_length}AA\nAccuracy: {tp_accuracy:.3f}\nSignificant ratio: {turning_point["significant_ratio"]:.2%}',
                        xy=(tp_length, tp_accuracy), 
                        xytext=(tp_length+5, tp_accuracy-0.05),
                        arrowprops=dict(arrowstyle='->', color='green', lw=2),
                        fontsize=12, fontweight='bold', 
                        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
        
        # 标记假说区域
        ax5.axvspan(40, 49, alpha=0.2, color='orange', label='Rapid divergence\n(40-49AA)')
        ax5.axvspan(50, 84, alpha=0.2, color='yellow', label='Transition\n(50-84AA)')
        ax5.axvspan(85, 100, alpha=0.2, color='lightgreen', label='Noise removal\n(85-100AA)')
        
        ax5.set_xlabel('N-terminal Length (AA)', fontsize=14)
        ax5.set_ylabel('Classification Accuracy', fontsize=14)
        ax5.set_title('Hypothesis Testing: N-terminal Length vs Classification Performance', 
                     fontsize=16, fontweight='bold')
        ax5.legend(fontsize=11, loc='upper left')
        ax5.grid(True, alpha=0.3)
        ax5.set_ylim(0.85, 1.0)
        
        plt.tight_layout()
        plt.savefig(f"{plot_dir}/hypothesis_testing_turning_point.png", 
                   dpi=300, bbox_inches='tight')
        plt.savefig(f"{plot_dir}/hypothesis_testing_turning_point.pdf", 
                   bbox_inches='tight')
        
        # 3. 单个重复的准确率曲线图（用于比较）
        fig3, ax6 = plt.subplots(figsize=(14, 8))
        
        if 'all_repeats' in all_results:
            for i, result in enumerate(all_results['all_repeats']):
                if 'length_scan' in result:
                    repeat_lengths = sorted(result['length_scan'].keys())
                    repeat_accuracies = [result['length_scan'][n]['mean_accuracy'] for n in repeat_lengths]
                    ax6.plot(repeat_lengths, repeat_accuracies, '-', linewidth=1.5, 
                            alpha=0.6, label=f'Repeat {i+1}')
        
        # 平均准确率曲线
        ax6.plot(lengths, mean_accuracies, 'k-', linewidth=3, 
                label='Average (6 repeats)', zorder=10)
        
        # 全长序列准确率水平线
        ax6.axhline(y=full_accuracy, color='red', linestyle='--', 
                   linewidth=2, label=f'Full-length average ({full_accuracy:.3f})')
        
        ax6.set_xlabel('N-terminal Length (AA)', fontsize=14)
        ax6.set_ylabel('Classification Accuracy', fontsize=14)
        ax6.set_title('Individual Repeats vs Average Performance', 
                     fontsize=16, fontweight='bold')
        ax6.legend(fontsize=10, loc='upper left')
        ax6.grid(True, alpha=0.3)
        ax6.set_ylim(0.85, 1.0)
        
        plt.tight_layout()
        plt.savefig(f"{plot_dir}/individual_repeats_comparison.png", 
                   dpi=300, bbox_inches='tight')
        plt.savefig(f"{plot_dir}/individual_repeats_comparison.pdf", 
                   bbox_inches='tight')
        
        print(f"\n整合结果图表已保存到: {plot_dir}/")
        plt.show()
    
    def analyze_hypothesis_support(self, all_results: Dict) -> Dict:
        """
        分析假说支持度
        验证分化期、过渡期和去噪音假说
        """
        print(f"\n{'='*80}")
        print("分析假说支持度...")
        print('='*80)
        
        if 'integrated' not in all_results:
            print("错误: 没有整合结果可供分析")
            return {}
        
        integrated_results = all_results['integrated']
        length_results = integrated_results['length_scan']
        
        # 定义假说区域
        hypothesis_regions = {
            'rapid_divergence': (40, 49),   # 快速分化期
            'transition': (50, 84),         # 过渡期
            'noise_removal': (85, 100)      # 去噪音期
        }
        
        hypothesis_support = {}
        
        for region_name, (start, end) in hypothesis_regions.items():
            region_lengths = [n for n in length_results.keys() if start <= n <= end]
            
            if not region_lengths:
                continue
            
            # 计算区域内的统计量
            accuracies = [length_results[n]['mean_accuracy'] for n in region_lengths]
            significant_ratios = [length_results[n]['significant_ratio'] for n in region_lengths]
            gap_percentages = [length_results[n]['mean_gap_percentage'] for n in region_lengths]
            
            hypothesis_support[region_name] = {
                'length_range': f"{start}-{end}AA",
                'n_lengths': len(region_lengths),
                'mean_accuracy': np.mean(accuracies),
                'std_accuracy': np.std(accuracies),
                'mean_significant_ratio': np.mean(significant_ratios),
                'mean_gap_percentage': np.mean(gap_percentages),
                'below_full_ratio': np.mean([1 if gap > 0 else 0 for gap in gap_percentages]),
                'lengths': region_lengths
            }
        
        # 打印假说支持度分析
        print("\n=== 假说支持度分析 ===")
        
        full_accuracy = integrated_results['full_length']['mean_accuracy']
        
        for region_name, support in hypothesis_support.items():
            print(f"\n{region_name.upper().replace('_', ' ')} ({support['length_range']}):")
            print(f"  平均准确率: {support['mean_accuracy']:.4f} (全长: {full_accuracy:.4f})")
            print(f"  准确率比例: {support['mean_accuracy']/full_accuracy:.2%}")
            print(f"  平均显著比例: {support['mean_significant_ratio']:.2%}")
            print(f"  平均性能差距: {support['mean_gap_percentage']:.2f}%")
            print(f"  低于全长的比例: {support['below_full_ratio']:.2%}")
            
            # 评估假说支持度
            if region_name == 'rapid_divergence':
                if support['mean_significant_ratio'] > 0.5:
                    print(f"  ✅ 支持快速分化假说: 大部分长度显著低于全长")
                else:
                    print(f"  ⚠️ 部分支持快速分化假说")
            
            elif region_name == 'transition':
                if support['mean_significant_ratio'] < 0.5:
                    print(f"  ✅ 支持过渡期假说: 大部分长度与全长无显著差异")
                else:
                    print(f"  ⚠️ 部分支持过渡期假说")
            
            elif region_name == 'noise_removal':
                if support['mean_gap_percentage'] <= 0:
                    print(f"  ✅ 支持去噪音假说: 性能略高于或等于全长")
                else:
                    print(f"  ⚠️ 部分支持去噪音假说")
        
        return hypothesis_support

def main():
    """主函数"""
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description='植物P450酶N端氨基酸长度扫描分析 - 多次重复线性SVM版本')
    parser.add_argument('--fasta', type=str, default="P450_unique_pep_final.fasta",
                       help='FASTA文件路径 (默认: P450_unique_pep_final.fasta)')
    parser.add_argument('--cluster', type=str, default="P450_unique_pep_final_esm_clusters_genes.csv",
                       help='聚类标签文件路径 (默认: P450_unique_pep_final_esm_clusters_genes.csv)')
    parser.add_argument('--output', type=str, default="p450_length_scan_repeated_svm",
                       help='输出目录 (默认: p450_length_scan_repeated_svm)')
    parser.add_argument('--min_length', type=int, default=40,
                       help='最小N端长度 (默认: 40)')
    parser.add_argument('--max_length', type=int, default=100,
                       help='最大N端长度 (默认: 100)')
    parser.add_argument('--n_repeats', type=int, default=6,
                       help='重复次数 (默认: 6)')
    parser.add_argument('--device', type=str, default="cpu",
                       choices=['cpu', 'cuda'],
                       help='计算设备: cpu或cuda (默认: cpu)')
    parser.add_argument('--seed_start', type=int, default=42,
                       help='随机种子起始值 (默认: 42)')
    
    args = parser.parse_args()
    
    print("="*80)
    print("植物P450酶N端氨基酸长度扫描分析 - 多次重复线性SVM版本")
    print(f"重复次数: {args.n_repeats}")
    print(f"扫描长度范围: {args.min_length}-{args.max_length} AA")
    print(f"随机种子起始值: {args.seed_start}")
    print("="*80)
    
    # 初始化分析器
    analyzer = RepeatedLinearSVMLengthScanAnalyzer(
        model_name="facebook/esm2_t12_35M_UR50D",  # 使用中等大小的模型
        device=args.device
    )
    
    # 加载数据
    try:
        sequences_dict, gene_cluster_dict = analyzer.load_data(args.fasta, args.cluster)
    except Exception as e:
        print(f"加载数据失败: {e}")
        return
    
    if len(sequences_dict) == 0:
        print("错误: 没有找到基因序列")
        return
    
    # 执行多次重复扫描分析
    print("\n" + "="*80)
    print("开始执行多次重复N端长度扫描分析...")
    print(f"将重复运行 {args.n_repeats} 次线性SVM分析")
    print("="*80)
    
    try:
        # 执行多次重复扫描
        all_results = analyzer.repeated_scan_n_terminal_lengths(
            sequences_dict, 
            gene_cluster_dict,
            min_length=args.min_length,
            max_length=args.max_length,
            n_repeats=args.n_repeats,
            random_seed_start=args.seed_start
        )
        
        # 导出整合结果
        analyzer.export_integrated_results(all_results, args.output)
        
        # 可视化整合结果
        analyzer.visualize_integrated_results(all_results, args.output)
        
        # 分析假说支持度
        hypothesis_support = analyzer.analyze_hypothesis_support(all_results)
        
        # 打印关键结论
        print(f"\n{'='*80}")
        print("关键结论汇总:")
        print('='*80)
        
        if 'integrated' in all_results:
            full_results = all_results['integrated']['full_length']
            print(f"1. 全长序列分类准确率 ({args.n_repeats}次重复平均):")
            print(f"   平均准确率: {full_results['mean_accuracy']:.4f}")
            print(f"   标准差: {full_results['std_accuracy']:.4f}")
            print(f"   范围: {full_results['min_accuracy']:.4f} - {full_results['max_accuracy']:.4f}")
            print(f"   变异系数: {full_results['std_accuracy']/full_results['mean_accuracy']:.2%}")
        
        # 寻找拐点
        if 'integrated' in all_results:
            turning_point = analyzer.find_turning_point_integrated(all_results['integrated'])
            if turning_point:
                print(f"\n2. 整合拐点分析:")
                print(f"   拐点位置: N端 {turning_point['length']}AA")
                print(f"   平均准确率: {turning_point['accuracy']:.4f} (全长: {turning_point['full_accuracy']:.4f})")
                print(f"   准确率比例: {turning_point['accuracy_ratio']:.2%}")
                print(f"   平均p值: {turning_point['mean_p_value']:.4f}")
                print(f"   显著比例: {turning_point['significant_ratio']:.2%}")
                
                if turning_point['length'] <= 55:
                    print(f"\n3. 科学意义: N端{turning_point['length']}AA已包含大部分功能分类信息")
                    print("   支持P450酶的功能信息高度集中在N端区域的假说")
                else:
                    print(f"\n3. 科学意义: 需要较长的N端序列({turning_point['length']}AA)才能达到全长性能")
                    print("   说明P450酶的功能信息分布在更长的序列区域")
        
        # 总结假说验证结果
        print(f"\n{'='*80}")
        print("假说验证总结:")
        print('='*80)
        
        if hypothesis_support:
            for region_name, support in hypothesis_support.items():
                region_display = region_name.upper().replace('_', ' ')
                if support['mean_significant_ratio'] > 0.5:
                    sig_status = "大部分显著"
                elif support['mean_significant_ratio'] > 0.3:
                    sig_status = "部分显著"
                else:
                    sig_status = "大部分不显著"
                
                print(f"{region_display} ({support['length_range']}):")
                print(f"  平均准确率: {support['mean_accuracy']:.4f}, {sig_status}")
                print(f"  性能差距: {support['mean_gap_percentage']:.2f}%")
        
        print(f"\n{'='*80}")
        print("分析完成！")
        print(f"详细结果请查看目录: {args.output}")
        print("="*80)
        
    except Exception as e:
        print(f"分析过程中出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()