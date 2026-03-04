#!/usr/bin/env python3
"""
增强版P450位置重要性分析脚本（修复JSON序列化问题）
"""

import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

from transformers import AutoTokenizer, AutoModel
from sklearn.svm import LinearSVC
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold
from scipy import stats
import scipy.stats as st
from tqdm import tqdm
import pickle
import json

class NumpyEncoder(json.JSONEncoder):
    """自定义JSON编码器，用于处理NumPy数据类型"""
    def default(self, obj):
        if isinstance(obj, (np.integer, np.int8, np.int16, np.int32, np.int64,
                          np.uint8, np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif pd.isna(obj):  # 处理NaN值
            return None
        else:
            return super().default(obj)

def convert_to_serializable(obj):
    """递归转换对象为JSON可序列化的类型"""
    if isinstance(obj, dict):
        return {key: convert_to_serializable(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(item) for item in obj]
    elif isinstance(obj, tuple):
        return tuple(convert_to_serializable(item) for item in obj)
    elif isinstance(obj, (np.integer, np.int8, np.int16, np.int32, np.int64,
                         np.uint8, np.uint16, np.uint32, np.uint64)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float16, np.float32, np.float64)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif pd.isna(obj):
        return None
    elif hasattr(obj, '__dict__'):  # 处理对象
        return convert_to_serializable(obj.__dict__)
    else:
        return obj

class EnhancedPositionImportanceAnalyzer:
    """增强版位置重要性分析器（6次重复+正选择整合）"""
    
    def __init__(self, model_name: str = "facebook/esm2_t33_650M_UR50D", 
                 device: str = "cpu"):
        self.device = device
        print(f"Loading ESM2 model: {model_name}...")
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModel.from_pretrained(model_name)
        self.model = self.model.to(self.device).eval()
        self.embedding_dim = self.model.embeddings.word_embeddings.embedding_dim
    
    def extract_embeddings(self, sequences: List[str], method: str = 'mean') -> np.ndarray:
        """提取序列嵌入向量"""
        all_embeddings = []
        batch_size = 2
        
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
                    
                    if method == 'mean':
                        seq_embeddings = last_hidden_state[j, 1:seq_length-1]
                        if len(seq_embeddings) > 0:
                            embedding = seq_embeddings.mean(dim=0).cpu().numpy()
                        else:
                            embedding = last_hidden_state[j, :seq_length].mean(dim=0).cpu().numpy()
                    elif method == 'cls':
                        embedding = last_hidden_state[j, 0].cpu().numpy()
                    else:
                        seq_embeddings = last_hidden_state[j, 1:seq_length-1]
                        if len(seq_embeddings) > 0:
                            embedding = seq_embeddings.mean(dim=0).cpu().numpy()
                        else:
                            embedding = last_hidden_state[j, :seq_length].mean(dim=0).cpu().numpy()
                    
                    all_embeddings.append(embedding)
        
        return np.array(all_embeddings)
    
    def extract_positional_embeddings(self, sequences: List[str], n_aa: int = 51) -> List[np.ndarray]:
        """提取每个氨基酸位置的独立嵌入向量"""
        print(f"Extracting positional embeddings for {len(sequences)} sequences...")
        
        all_positional_embeddings = []
        batch_size = 2
        
        with torch.no_grad():
            for i in tqdm(range(0, len(sequences), batch_size), desc="Processing sequences"):
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
                
                for j, seq in enumerate(batch_sequences):
                    seq_length = len(seq) + 2  # Include [CLS] and [SEP]
                    seq_embeddings = last_hidden_state[j, 1:seq_length-1].cpu().numpy()
                    
                    # Limit to first n_aa positions
                    if len(seq_embeddings) > n_aa:
                        seq_embeddings = seq_embeddings[:n_aa]
                    all_positional_embeddings.append(seq_embeddings)
        
        return all_positional_embeddings
    
    def single_repeat_analysis(self, sequences: List[str], labels: List[str], 
                              n_aa: int = 51, random_state: int = 42) -> Dict:
        """单次重复的位置重要性分析"""
        print(f"Repeat analysis (Random seed: {random_state})...")
        
        # 设置随机种子
        np.random.seed(random_state)
        
        # 1. 准备数据
        X_mean = self.extract_embeddings(sequences, method='mean')
        y = LabelEncoder().fit_transform(labels)
        
        # 2. 数据预处理
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_mean)
        
        # PCA降维（提高稳定性）
        n_components = min(50, X_scaled.shape[0] - 1, X_scaled.shape[1])
        pca = PCA(n_components=n_components, random_state=random_state)
        X_pca = pca.fit_transform(X_scaled)
        
        # 3. 训练线性SVM（5折交叉验证）
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)
        fold_accuracies = []
        fold_weights = []
        
        for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X_pca, y)):
            X_train, X_test = X_pca[train_idx], X_pca[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            
            # 训练SVM
            svm = LinearSVC(
                C=1.0,
                class_weight='balanced',
                max_iter=10000,
                random_state=random_state,
                dual=False,
                tol=1e-4
            )
            svm.fit(X_train, y_train)
            
            # 计算准确率
            accuracy = svm.score(X_test, y_test)
            fold_accuracies.append(accuracy)
            
            # 保存权重
            fold_weights.append(svm.coef_[0].copy())
        
        # 4. 计算位置重要性
        positional_embeddings = self.extract_positional_embeddings(sequences, n_aa)
        position_importance = self.calculate_position_importance(
            fold_weights, positional_embeddings, pca, scaler, n_aa
        )
        
        # 5. 计算稳定性指标
        stability_scores = self.calculate_stability_metrics(fold_weights, pca, scaler, 
                                                          positional_embeddings, n_aa)
        
        return {
            'mean_accuracy': float(np.mean(fold_accuracies)),
            'std_accuracy': float(np.std(fold_accuracies)),
            'fold_accuracies': [float(acc) for acc in fold_accuracies],
            'fold_weights': [w.tolist() for w in fold_weights],
            'position_importance': position_importance,
            'stability_scores': stability_scores,
            'random_state': random_state,
            'pca_explained_variance': float(pca.explained_variance_ratio_.sum()),
            'n_pca_components': n_components
        }
    
    def calculate_position_importance(self, fold_weights: List[np.ndarray], 
                                    positional_embeddings: List[np.ndarray],
                                    pca: PCA, scaler: StandardScaler, 
                                    n_aa: int) -> Dict:
        """基于SVM权重计算位置重要性"""
        print("Calculating position importance...")
        
        # 1. 计算每个折叠的原始空间权重
        original_weights_list = []
        for w_pca in fold_weights:
            # 转换回原始空间: w_original = pca.components_.T @ w_pca
            w_original = pca.components_.T @ w_pca
            original_weights_list.append(w_original)
        
        # 2. 计算每个位置的重要性
        n_positions = n_aa
        position_scores_all_folds = np.zeros((len(fold_weights), n_positions))
        
        for fold_idx, w_original in enumerate(original_weights_list):
            for pos in range(n_positions):
                pos_projections = []
                
                for seq_idx, seq_embeds in enumerate(positional_embeddings):
                    if pos < len(seq_embeds):
                        pos_embed = seq_embeds[pos].flatten()
                        # 标准化（与训练时一致）
                        if len(pos_embed) == w_original.shape[0]:
                            pos_embed_scaled = scaler.transform(pos_embed.reshape(1, -1))
                            projection = np.abs(np.dot(w_original, pos_embed_scaled.flatten()))
                            pos_projections.append(projection)
                
                if pos_projections:
                    position_scores_all_folds[fold_idx, pos] = np.mean(pos_projections)
        
        # 3. 计算统计量
        position_importance = {
            'mean_scores': position_scores_all_folds.mean(axis=0).tolist(),
            'std_scores': position_scores_all_folds.std(axis=0).tolist(),
            'fold_scores': position_scores_all_folds.tolist(),
            'coefficient_of_variation': (position_scores_all_folds.std(axis=0) / 
                                       (position_scores_all_folds.mean(axis=0) + 1e-10)).tolist()
        }
        
        # 4. 标准化重要性分数（0-1范围）
        mean_scores = np.array(position_importance['mean_scores'])
        if mean_scores.max() > 0:
            position_importance['normalized_scores'] = (mean_scores / mean_scores.max()).tolist()
        else:
            position_importance['normalized_scores'] = mean_scores.tolist()
        
        return position_importance
    
    def calculate_stability_metrics(self, fold_weights: List[np.ndarray], 
                                  pca: PCA, scaler: StandardScaler,
                                  positional_embeddings: List[np.ndarray], 
                                  n_aa: int) -> Dict:
        """计算权重稳定性指标"""
        n_folds = len(fold_weights)
        
        # 1. 计算权重之间的相关性
        weight_correlations = []
        for i in range(n_folds):
            for j in range(i+1, n_folds):
                corr = np.corrcoef(fold_weights[i], fold_weights[j])[0, 1]
                weight_correlations.append(float(corr))
        
        # 2. 计算位置重要性排名的一致性
        position_rankings = []
        for w_pca in fold_weights:
            w_original = pca.components_.T @ w_pca
            pos_scores = np.zeros(n_aa)
            
            for pos in range(n_aa):
                pos_projections = []
                for seq_embeds in positional_embeddings:
                    if pos < len(seq_embeds):
                        pos_embed = seq_embeds[pos].flatten()
                        if len(pos_embed) == w_original.shape[0]:
                            pos_embed_scaled = scaler.transform(pos_embed.reshape(1, -1))
                            projection = np.abs(np.dot(w_original, pos_embed_scaled.flatten()))
                            pos_projections.append(projection)
                
                if pos_projections:
                    pos_scores[pos] = np.mean(pos_projections)
            
            # 获取位置排名
            ranking = np.argsort(-pos_scores)  # 降序排列
            position_rankings.append(ranking)
        
        # 3. 计算排名一致性（Kendall's W）
        if len(position_rankings) > 1:
            # 转换为排名矩阵
            rank_matrix = np.array(position_rankings)
            
            # 计算Kendall's W一致性系数
            m, n = rank_matrix.shape
            sum_ranks = rank_matrix.sum(axis=0)
            mean_ranks = sum_ranks / m
            
            # 计算每个位置的排名平方差
            S = np.sum((sum_ranks - mean_ranks.sum()/n) ** 2)
            
            # 计算W
            W = 12 * S / (m**2 * (n**3 - n))
        else:
            W = 1.0
        
        return {
            'mean_weight_correlation': float(np.mean(weight_correlations) if weight_correlations else 0),
            'std_weight_correlation': float(np.std(weight_correlations) if weight_correlations else 0),
            'kendalls_W': float(W),
            'n_folds': n_folds
        }
    
    def six_repeats_analysis(self, sequences: List[str], labels: List[str], 
                            n_aa: int = 51) -> Dict:
        """运行6次重复实验并整合结果"""
        print("="*80)
        print("Running 6-Repeat Position Importance Analysis")
        print("="*80)
        
        all_repeats = []
        random_seeds = [42, 123, 456, 789, 101112, 131415]  # 6个不同的随机种子
        
        for i, seed in enumerate(random_seeds):
            print(f"\nRepeat {i+1}/6 (Random seed: {seed})")
            
            try:
                repeat_result = self.single_repeat_analysis(
                    sequences, labels, n_aa, random_state=seed
                )
                all_repeats.append(repeat_result)
                
                # 保存单次结果
                self.save_single_repeat_result(repeat_result, i+1, seed)
                
            except Exception as e:
                print(f"Repeat {i+1} failed: {e}")
                continue
        
        # 整合结果
        integrated_results = self.integrate_repeat_results(all_repeats, n_aa)
        
        return {
            'all_repeats': all_repeats,
            'integrated': integrated_results,
            'n_successful_repeats': len(all_repeats),
            'random_seeds_used': random_seeds[:len(all_repeats)]
        }
    
    def integrate_repeat_results(self, all_repeats: List[Dict], n_aa: int) -> Dict:
        """整合6次重复的结果"""
        print(f"\nIntegrating {len(all_repeats)} repeat results...")
        
        if not all_repeats:
            raise ValueError("No repeat results to integrate")
        
        n_repeats = len(all_repeats)
        
        # 1. 整合位置重要性分数
        position_scores_matrix = np.zeros((n_repeats, n_aa))
        accuracies = []
        
        for i, repeat in enumerate(all_repeats):
            mean_scores = repeat['position_importance']['mean_scores']
            if len(mean_scores) >= n_aa:
                position_scores_matrix[i, :] = mean_scores[:n_aa]
            else:
                position_scores_matrix[i, :len(mean_scores)] = mean_scores
            accuracies.append(repeat['mean_accuracy'])
        
        # 2. 计算整合统计量
        integrated_position_importance = {
            'mean_across_repeats': position_scores_matrix.mean(axis=0).tolist(),
            'std_across_repeats': position_scores_matrix.std(axis=0).tolist(),
            'median_across_repeats': np.median(position_scores_matrix, axis=0).tolist(),
            'min_across_repeats': position_scores_matrix.min(axis=0).tolist(),
            'max_across_repeats': position_scores_matrix.max(axis=0).tolist(),
            'coefficient_of_variation': (position_scores_matrix.std(axis=0) / 
                                       (position_scores_matrix.mean(axis=0) + 1e-10)).tolist()
        }
        
        # 计算重复之间的相关性
        if n_repeats > 1:
            corr_matrix = np.corrcoef(position_scores_matrix)
            integrated_position_importance['repeat_correlations'] = corr_matrix[np.triu_indices(n_repeats, k=1)].tolist()
        else:
            integrated_position_importance['repeat_correlations'] = []
        
        # 3. 识别稳定关键位置（在多次重复中一致重要）
        # 计算每个位置的排名稳定性
        position_rankings = np.argsort(-position_scores_matrix, axis=1)  # 每行（每个重复）的位置排名
        
        # 计算每个位置的平均排名
        mean_ranks = np.zeros(n_aa)
        for pos in range(n_aa):
            ranks = []
            for repeat_idx in range(n_repeats):
                rank = np.where(position_rankings[repeat_idx] == pos)[0][0]
                ranks.append(rank)
            mean_ranks[pos] = np.mean(ranks)
        
        # 4. 关键位置标准：平均排名在前20%且变异系数小于0.5
        top_percentage = 0.2
        n_top = int(n_aa * top_percentage)
        cv_threshold = 0.5
        
        # 按平均排名排序
        sorted_positions = np.argsort(mean_ranks)
        candidate_positions = sorted_positions[:n_top]
        
        # 筛选变异系数低的
        key_positions = []
        for pos in candidate_positions:
            cv = integrated_position_importance['coefficient_of_variation'][pos]
            if cv < cv_threshold:
                key_positions.append(int(pos + 1))  # 转换为1-based索引
        
        integrated_position_importance['key_positions'] = key_positions
        integrated_position_importance['mean_ranks'] = mean_ranks.tolist()
        
        # 5. 整合准确率
        integrated_accuracy = {
            'mean_accuracy': float(np.mean(accuracies)),
            'std_accuracy': float(np.std(accuracies)),
            'min_accuracy': float(np.min(accuracies)),
            'max_accuracy': float(np.max(accuracies)),
            'accuracy_values': [float(acc) for acc in accuracies]
        }
        
        # 6. 整合稳定性指标
        stability_metrics = {
            'mean_weight_correlation': float(np.mean([r['stability_scores']['mean_weight_correlation'] 
                                               for r in all_repeats])),
            'mean_kendalls_W': float(np.mean([r['stability_scores']['kendalls_W'] 
                                       for r in all_repeats])),
            'repeat_consistency': float(np.mean(integrated_position_importance['repeat_correlations'])) if integrated_position_importance['repeat_correlations'] else 0.0
        }
        
        return {
            'position_importance': integrated_position_importance,
            'accuracy': integrated_accuracy,
            'stability': stability_metrics,
            'n_repeats': n_repeats
        }
    
    def compare_with_positive_selection(self, key_positions: List[int], 
                                       positive_selection_positions: List[int],
                                       n_aa: int = 51) -> Dict:
        """
        比较LinearSVM关键位置与正选择位点
        
        Args:
            key_positions: LinearSVM识别的关键位置（1-based）
            positive_selection_positions: 正选择位点（1-based）
            n_aa: 分析的总位置数
        """
        print(f"\nComparing with positive selection positions...")
        print(f"Positive selection positions: {positive_selection_positions}")
        print(f"LinearSVM key positions: {key_positions}")
        
        # 转换为集合以便操作
        svm_set = set(key_positions)
        ps_set = set(positive_selection_positions)
        
        # 1. 计算重叠
        overlap = svm_set.intersection(ps_set)
        n_overlap = len(overlap)
        
        # 2. 统计检验（超几何检验）
        # 原假设：重叠是随机的
        # 备择假设：重叠显著多于随机
        total_positions = n_aa
        n_svm = len(svm_set)
        n_ps = len(ps_set)
        
        # 超几何分布p值
        from scipy.stats import hypergeom
        # 计算观察到n_overlap或更多重叠的概率
        p_value = 1 - hypergeom.cdf(n_overlap - 1, total_positions, n_ps, n_svm)
        
        # 3. 计算其他统计量
        # 重叠比例
        overlap_ratio_svm = n_overlap / n_svm if n_svm > 0 else 0
        overlap_ratio_ps = n_overlap / n_ps if n_ps > 0 else 0
        
        # Jaccard相似系数
        union_size = len(svm_set.union(ps_set))
        jaccard = n_overlap / union_size if union_size > 0 else 0
        
        # 4. 计算每个位置类别的统计
        categories = {
            'both_important': sorted(list(overlap)),
            'svm_only': sorted(list(svm_set - ps_set)),
            'ps_only': sorted(list(ps_set - svm_set)),
            'neither': sorted(list(set(range(1, n_aa+1)) - svm_set - ps_set))
        }
        
        comparison_results = {
            'overlap_positions': sorted(list(overlap)),
            'n_overlap': int(n_overlap),
            'p_value': float(p_value),
            'overlap_ratio_svm': float(overlap_ratio_svm),
            'overlap_ratio_ps': float(overlap_ratio_ps),
            'jaccard_similarity': float(jaccard),
            'hypergeometric_test': {
                'total_positions': int(total_positions),
                'n_svm_positions': int(n_svm),
                'n_ps_positions': int(n_ps),
                'expected_overlap': float(n_svm * n_ps / total_positions),
                'observed_overlap': int(n_overlap),
                'p_value': float(p_value)
            },
            'categories': categories,
            'summary': {
                'significance': 'Significant' if p_value < 0.05 else 'Not significant',
                'interpretation': self.interpret_overlap_results(p_value, overlap_ratio_svm)
            }
        }
        
        # 打印结果
        print(f"\n=== Comparison Results ===")
        print(f"Overlap positions: {sorted(list(overlap))}")
        print(f"Number of overlapping positions: {n_overlap}")
        print(f"Hypergeometric test p-value: {p_value:.6f}")
        print(f"Significance: {'Significant (p < 0.05)' if p_value < 0.05 else 'Not significant'}")
        print(f"Overlap ratio in SVM positions: {overlap_ratio_svm:.2%}")
        print(f"Overlap ratio in PS positions: {overlap_ratio_ps:.2%}")
        print(f"Jaccard similarity: {jaccard:.4f}")
        
        return comparison_results
    
    def interpret_overlap_results(self, p_value: float, overlap_ratio: float) -> str:
        """解释重叠结果的生物学意义"""
        interpretations = []
        
        if p_value < 0.05:
            interpretations.append("SVM关键位置与正选择位点的重叠具有统计学显著性")
            
            if overlap_ratio > 0.3:
                interpretations.append("较高的重叠比例提示这些位置可能在功能分类和适应性进化中都起关键作用")
            else:
                interpretations.append("有限的重叠提示功能分类和适应性进化可能涉及不同的残基")
        else:
            interpretations.append("SVM关键位置与正选择位点的重叠不显著")
            interpretations.append("功能分类的关键残基可能与适应性进化的残基不同")
        
        if p_value < 0.01 and overlap_ratio > 0.4:
            interpretations.append("强烈支持正选择位点在功能分类中起关键作用的假说")
        
        return "; ".join(interpretations)
    
    def save_single_repeat_result(self, repeat_result: Dict, repeat_idx: int, random_seed: int):
        """保存单次重复结果"""
        output_dir = f"repeat_{repeat_idx}_seed_{random_seed}"
        Path(output_dir).mkdir(exist_ok=True)
        
        # 保存位置重要性数据
        importance_df = pd.DataFrame({
            'position': range(1, len(repeat_result['position_importance']['mean_scores']) + 1),
            'importance_score': repeat_result['position_importance']['mean_scores'],
            'std_score': repeat_result['position_importance']['std_scores'],
            'normalized_score': repeat_result['position_importance']['normalized_scores']
        })
        importance_df.to_csv(f"{output_dir}/position_importance.csv", index=False)
        
        # 保存元数据
        metadata = {
            'repeat_index': repeat_idx,
            'random_seed': random_seed,
            'mean_accuracy': repeat_result['mean_accuracy'],
            'std_accuracy': repeat_result['std_accuracy'],
            'pca_explained_variance': repeat_result['pca_explained_variance'],
            'n_pca_components': repeat_result['n_pca_components'],
            'stability_metrics': repeat_result['stability_scores']
        }
        
        # 转换元数据为可序列化格式
        metadata_serializable = convert_to_serializable(metadata)
        
        with open(f"{output_dir}/metadata.json", 'w') as f:
            json.dump(metadata_serializable, f, indent=2, cls=NumpyEncoder)
    
    def visualize_integrated_results(self, integrated_results: Dict, 
                                   positive_selection_positions: List[int] = None,
                                   output_dir: str = "integrated_results"):
        """可视化整合结果"""
        print("\nGenerating integrated visualizations...")
        Path(output_dir).mkdir(exist_ok=True)
        plot_dir = Path(output_dir) / "plots"
        plot_dir.mkdir(exist_ok=True)
        
        plt.style.use('seaborn-v0_8-whitegrid')
        sns.set_palette("husl")
        
        # 获取数据
        pos_importance = integrated_results['position_importance']
        positions = np.arange(1, len(pos_importance['mean_across_repeats']) + 1)
        mean_scores = np.array(pos_importance['mean_across_repeats'])
        std_scores = np.array(pos_importance['std_across_repeats'])
        key_positions = pos_importance.get('key_positions', [])
        
        # 1. 主要位置重要性图（带误差条）
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
        
        # 图1: 位置重要性条形图（带标准差）
        bars = ax1.bar(positions, mean_scores, yerr=std_scores, 
                      capsize=3, alpha=0.7, color='steelblue', error_kw={'elinewidth': 1})
        ax1.set_xlabel('Amino Acid Position', fontsize=12)
        ax1.set_ylabel('Importance Score', fontsize=12)
        ax1.set_title('Position Importance (6 Repeats Average)', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='y')
        
        # 标记关键位置
        for pos in key_positions:
            idx = pos - 1
            ax1.annotate(f"Pos{pos}", xy=(pos, mean_scores[idx]),
                        xytext=(0, 10), textcoords="offset points",
                        ha='center', fontsize=9, fontweight='bold', color='red')
        
        # 标记正选择位点（如果有）
        if positive_selection_positions:
            for pos in positive_selection_positions:
                if 1 <= pos <= len(mean_scores):
                    idx = pos - 1
                    ax1.annotate(f"PS{pos}", xy=(pos, mean_scores[idx]),
                                xytext=(0, -15), textcoords="offset points",
                                ha='center', fontsize=8, fontweight='bold', 
                                color='green', bbox=dict(boxstyle='round,pad=0.2', 
                                                       facecolor='lightgreen', alpha=0.7))
        
        ax1.legend(['Importance ± SD', 'SVM Key Position', 'Positive Selection'], fontsize=10)
        
        # 图2: 位置重要性热图（跨重复）
        if 'all_repeats' in integrated_results:
            n_repeats = len(integrated_results['all_repeats'])
            heatmap_data = np.zeros((n_repeats, len(positions)))
            
            for i, repeat in enumerate(integrated_results['all_repeats'][:n_repeats]):
                scores = repeat['position_importance']['mean_scores']
                heatmap_data[i, :len(scores)] = scores[:len(positions)]
            
            im = ax2.imshow(heatmap_data, aspect='auto', cmap='viridis',
                           extent=[0.5, len(positions)+0.5, n_repeats+0.5, 0.5])
            ax2.set_xlabel('Amino Acid Position', fontsize=12)
            ax2.set_ylabel('Repeat', fontsize=12)
            ax2.set_title('Position Importance Across 6 Repeats', fontsize=14, fontweight='bold')
            ax2.set_xticks(range(1, len(positions)+1, 5))
            ax2.set_yticks(range(1, n_repeats+1))
            
            plt.colorbar(im, ax=ax2, label='Importance Score')
        
        # 图3: 累积重要性分布
        sorted_indices = np.argsort(-mean_scores)
        cumulative_importance = np.cumsum(mean_scores[sorted_indices]) / np.sum(mean_scores)
        
        ax3.plot(range(1, len(positions)+1), cumulative_importance, 
                'b-', linewidth=2, marker='o', markersize=4)
        ax3.set_xlabel('Number of Positions (Sorted by Importance)', fontsize=12)
        ax3.set_ylabel('Cumulative Importance', fontsize=12)
        ax3.set_title('Cumulative Importance Distribution', fontsize=14, fontweight='bold')
        ax3.grid(True, alpha=0.3)
        
        # 标记关键阈值
        for threshold in [0.5, 0.8, 0.9]:
            idx = np.where(cumulative_importance >= threshold)[0]
            if len(idx) > 0:
                n_pos = idx[0] + 1
                ax3.axvline(x=n_pos, color='red' if threshold == 0.5 else 
                           'orange' if threshold == 0.8 else 'green',
                           linestyle='--', alpha=0.7)
                ax3.annotate(f'{threshold:.0%} at {n_pos} positions',
                            xy=(n_pos, threshold), xytext=(10, 10),
                            textcoords="offset points", fontsize=10,
                            color='red' if threshold == 0.5 else 
                            'orange' if threshold == 0.8 else 'green')
        
        # 图4: 准确率分布（箱线图）
        if 'accuracy' in integrated_results:
            acc_data = integrated_results['accuracy']['accuracy_values']
            
            # 箱线图
            bp = ax4.boxplot([acc_data], patch_artist=True)
            bp['boxes'][0].set_facecolor('lightblue')
            
            # 添加散点显示每个重复
            for i, acc in enumerate(acc_data):
                ax4.scatter(1, acc, color='red', s=50, zorder=3)
            
            ax4.set_ylabel('Accuracy', fontsize=12)
            ax4.set_title(f'Accuracy Distribution (Mean: {np.mean(acc_data):.3f} ± {np.std(acc_data):.3f})', 
                         fontsize=14, fontweight='bold')
            ax4.set_xticks([1])
            ax4.set_xticklabels(['6 Repeats'])
            ax4.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(f"{plot_dir}/integrated_position_importance.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{plot_dir}/integrated_position_importance.pdf", bbox_inches='tight')
        
        # 2. 与正选择位点的比较图（如果有）
        if positive_selection_positions:
            fig2, ax5 = plt.subplots(figsize=(14, 8))
            
            # 绘制重要性曲线
            ax5.plot(positions, mean_scores, 'b-', linewidth=2, label='LinearSVM Importance')
            ax5.fill_between(positions, mean_scores - std_scores, mean_scores + std_scores,
                           alpha=0.2, color='blue', label='± SD')
            
            # 标记正选择位点
            for pos in positive_selection_positions:
                if 1 <= pos <= len(mean_scores):
                    idx = pos - 1
                    ax5.scatter(pos, mean_scores[idx], color='green', s=100, 
                              zorder=5, label='Positive Selection' if pos == positive_selection_positions[0] else '')
                    ax5.annotate(f'PS{pos}', xy=(pos, mean_scores[idx]),
                                xytext=(0, 15), textcoords="offset points",
                                ha='center', fontsize=9, fontweight='bold', color='green')
            
            # 标记SVM关键位置
            for pos in key_positions:
                idx = pos - 1
                ax5.scatter(pos, mean_scores[idx], color='red', s=150, marker='*',
                          zorder=6, label='SVM Key Position' if pos == key_positions[0] else '')
                ax5.annotate(f'Key{pos}', xy=(pos, mean_scores[idx]),
                            xytext=(0, -20), textcoords="offset points",
                            ha='center', fontsize=9, fontweight='bold', color='red')
            
            ax5.set_xlabel('Amino Acid Position', fontsize=12)
            ax5.set_ylabel('Importance Score', fontsize=12)
            ax5.set_title('LinearSVM Importance vs Positive Selection Positions', 
                         fontsize=14, fontweight='bold')
            ax5.legend(fontsize=10, loc='upper right')
            ax5.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(f"{plot_dir}/svm_vs_positive_selection.png", dpi=300, bbox_inches='tight')
            plt.savefig(f"{plot_dir}/svm_vs_positive_selection.pdf", bbox_inches='tight')
        
        # 3. 稳定性分析图
        fig3, ax6 = plt.subplots(figsize=(12, 8))
        
        # 绘制变异系数
        cv = np.array(pos_importance['coefficient_of_variation'])
        ax6.bar(positions, cv, alpha=0.7, color='orange')
        ax6.axhline(y=0.5, color='red', linestyle='--', linewidth=2, label='CV Threshold (0.5)')
        
        # 标记稳定性好的位置（CV < 0.5）
        stable_positions = positions[cv < 0.5]
        for pos in stable_positions:
            ax6.annotate(f"{int(pos)}", xy=(pos, cv[int(pos)-1]),
                        xytext=(0, 5), textcoords="offset points",
                        ha='center', fontsize=8, color='darkgreen')
        
        ax6.set_xlabel('Amino Acid Position', fontsize=12)
        ax6.set_ylabel('Coefficient of Variation', fontsize=12)
        ax6.set_title('Position Importance Stability (Lower CV = More Stable)', 
                     fontsize=14, fontweight='bold')
        ax6.legend(fontsize=10)
        ax6.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(f"{plot_dir}/position_stability.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{plot_dir}/position_stability.pdf", bbox_inches='tight')
        
        print(f"Visualizations saved to: {plot_dir}/")
        plt.show()
    
    def generate_comprehensive_report(self, integrated_results: Dict, 
                                    comparison_results: Dict = None,
                                    output_file: str = "comprehensive_analysis_report.txt"):
        """生成综合分析报告"""
        print("Generating comprehensive analysis report...")
        
        # 确保输出目录存在
        output_dir = Path(output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("COMPREHENSIVE P450 N-TERMINAL POSITION IMPORTANCE ANALYSIS REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            # 1. 分析摘要
            f.write("1. EXECUTIVE SUMMARY\n")
            f.write("-" * 50 + "\n")
            f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Analysis Method: LinearSVM on ESM2 embeddings (6 repeats)\n")
            f.write(f"Amino Acid Region: Positions 1-{len(integrated_results['position_importance']['mean_across_repeats'])}\n")
            f.write(f"Successful Repeats: {integrated_results.get('n_repeats', 0)}/6\n")
            f.write(f"Mean Classification Accuracy: {integrated_results['accuracy']['mean_accuracy']:.4f} "
                   f"(±{integrated_results['accuracy']['std_accuracy']:.4f})\n")
            f.write(f"Accuracy Range: {integrated_results['accuracy']['min_accuracy']:.4f} - "
                   f"{integrated_results['accuracy']['max_accuracy']:.4f}\n")
            f.write(f"Coefficient of Variation (Accuracy): "
                   f"{integrated_results['accuracy']['std_accuracy']/integrated_results['accuracy']['mean_accuracy']:.2%}\n\n")
            
            # 2. 关键位置识别
            f.write("2. KEY POSITION IDENTIFICATION\n")
            f.write("-" * 50 + "\n")
            
            key_positions = integrated_results['position_importance'].get('key_positions', [])
            f.write(f"Number of Key Positions Identified: {len(key_positions)}\n")
            f.write(f"Key Positions (1-based): {sorted(key_positions)}\n\n")
            
            # 关键位置的统计特征
            if key_positions:
                f.write("Top 10 Most Important Positions:\n")
                mean_scores = np.array(integrated_results['position_importance']['mean_across_repeats'])
                sorted_indices = np.argsort(-mean_scores)[:10]  # 降序排列，取前10
                
                f.write("{:<8} {:<15} {:<15} {:<15}\n".format(
                    "Position", "Mean Score", "Std Score", "CV"
                ))
                f.write("-" * 60 + "\n")
                
                for idx in sorted_indices:
                    pos = idx + 1
                    mean_score = mean_scores[idx]
                    std_score = integrated_results['position_importance']['std_across_repeats'][idx]
                    cv = integrated_results['position_importance']['coefficient_of_variation'][idx]
                    
                    f.write("{:<8} {:<15.4f} {:<15.4f} {:<15.4f}\n".format(
                        f"Pos{pos}", mean_score, std_score, cv
                    ))
                f.write("\n")
            
            # 3. 稳定性分析
            f.write("3. ANALYSIS STABILITY AND ROBUSTNESS\n")
            f.write("-" * 50 + "\n")
            
            stability = integrated_results.get('stability', {})
            f.write(f"Mean Weight Correlation Across Folds: {stability.get('mean_weight_correlation', 0):.4f}\n")
            f.write(f"Mean Kendall's W (Ranking Consistency): {stability.get('mean_kendalls_W', 0):.4f}\n")
            f.write(f"Repeat-to-Repeat Consistency: {stability.get('repeat_consistency', 0):.4f}\n\n")
            
            # 4. 与正选择位点比较（如果有）
            if comparison_results:
                f.write("4. COMPARISON WITH POSITIVE SELECTION POSITIONS\n")
                f.write("-" * 50 + "\n")
                
                f.write(f"Positive Selection Positions Provided: {comparison_results['hypergeometric_test']['n_ps_positions']}\n")
                f.write(f"LinearSVM Key Positions: {comparison_results['hypergeometric_test']['n_svm_positions']}\n")
                f.write(f"Overlap Positions: {comparison_results['n_overlap']}\n")
                f.write(f"Overlap Positions List: {sorted(comparison_results['overlap_positions'])}\n")
                f.write(f"Jaccard Similarity Coefficient: {comparison_results['jaccard_similarity']:.4f}\n")
                f.write(f"Hypergeometric Test p-value: {comparison_results['p_value']:.6f}\n")
                f.write(f"Statistical Significance: {comparison_results['summary']['significance']}\n\n")
                
                # 重叠的生物学解释
                f.write("Biological Interpretation of Overlap:\n")
                f.write(f"{comparison_results['summary']['interpretation']}\n\n")
                
                # 分类统计
                f.write("Position Classification Summary:\n")
                categories = comparison_results['categories']
                f.write(f"  Positions Important for Both: {len(categories['both_important'])} positions\n")
                if categories['both_important']:
                    f.write(f"    List: {categories['both_important']}\n")
                
                f.write(f"  SVM-Only Positions: {len(categories['svm_only'])} positions\n")
                if categories['svm_only'] and len(categories['svm_only']) <= 20:
                    f.write(f"    List: {categories['svm_only']}\n")
                
                f.write(f"  Positive-Selection-Only Positions: {len(categories['ps_only'])} positions\n")
                if categories['ps_only'] and len(categories['ps_only']) <= 20:
                    f.write(f"    List: {categories['ps_only']}\n")
                
                f.write(f"  Positions Not Important for Either: {len(categories['neither'])} positions\n\n")
            
            # 5. 生物学意义
            f.write("5. BIOLOGICAL SIGNIFICANCE\n")
            f.write("-" * 50 + "\n")
            
            if key_positions:
                f.write("Functional Hypotheses for Key Positions:\n")
                
                # 根据位置分组
                early_positions = [p for p in key_positions if p <= 20]
                mid_positions = [p for p in key_positions if 21 <= p <= 40]
                late_positions = [p for p in key_positions if p > 40]
                
                if early_positions:
                    f.write(f"  Early N-terminal (1-20): {sorted(early_positions)}\n")
                    f.write("    Potential roles: Signal peptide, membrane targeting, initial folding\n\n")
                
                if mid_positions:
                    f.write(f"  Mid N-terminal (21-40): {sorted(mid_positions)}\n")
                    f.write("    Potential roles: Substrate channel entry, structural motifs\n\n")
                
                if late_positions:
                    f.write(f"  Late N-terminal (>40): {sorted(late_positions)}\n")
                    f.write("    Potential roles: Active site proximity, cofactor binding\n\n")
            
            # 6. 实验验证建议
            f.write("6. EXPERIMENTAL VALIDATION RECOMMENDATIONS\n")
            f.write("-" * 50 + "\n")
            
            if key_positions:
                top_5 = sorted(key_positions)[:5]
                f.write(f"Priority Targets for Site-Directed Mutagenesis:\n")
                f.write(f"  1. Top 5 Key Positions: {top_5}\n")
                f.write("  2. Recommended Experiments:\n")
                f.write("     a. Alanine scanning of key positions\n")
                f.write("     b. Truncation analysis at key positions\n")
                f.write("     c. Enzyme activity assays with mutants\n")
                f.write("     d. Structural analysis (if structures available)\n")
                f.write("     e. Subcellular localization studies\n\n")
            
            # 7. 方法学评估
            f.write("7. METHODOLOGICAL ASSESSMENT\n")
            f.write("-" * 50 + "\n")
            
            f.write("Strengths of the Analysis:\n")
            f.write("  1. 6-fold repetition ensures robustness\n")
            f.write("  2. LinearSVM provides interpretable weights\n")
            f.write("  3. ESM2 embeddings capture evolutionary and structural information\n")
            f.write("  4. Statistical validation through hypergeometric testing\n\n")
            
            f.write("Limitations and Considerations:\n")
            f.write("  1. Linear model may miss non-linear interactions\n")
            f.write("  2. Results are specific to the classification task\n")
            f.write("  3. Sample size may affect stability of estimates\n")
            f.write("  4. Position importance is relative, not absolute\n\n")
            
            # 8. 结论
            f.write("8. CONCLUSIONS\n")
            f.write("-" * 50 + "\n")
            
            f.write("This comprehensive analysis identified key amino acid positions in the P450 N-terminal region\n")
            f.write("that are important for functional classification. The 6-repeat design provides robust estimates\n")
            f.write("of position importance, and comparison with positive selection positions offers evolutionary context.\n")
            
            if comparison_results and comparison_results['p_value'] < 0.05:
                f.write("The significant overlap between LinearSVM key positions and positive selection positions\n")
                f.write("suggests that residues important for functional classification are also under positive selection,\n")
                f.write("supporting their biological importance in P450 evolution and function.\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write("END OF REPORT\n")
            f.write("=" * 80 + "\n")
        
        print(f"Comprehensive report generated: {output_file}")

def main():
    """主函数"""
    print("=" * 80)
    print("ENHANCED P450 POSITION IMPORTANCE ANALYSIS (6 REPEATS)")
    print("=" * 80)
    
    # 配置参数
    FASTA_FILE = "P450_unique_pep_final.fasta"
    CLUSTER_FILE = "P450_unique_pep_final_esm_clusters_genes.csv"
    OUTPUT_DIR = "enhanced_position_importance_6repeats"
    
    # 正选择位点（从您提供的列表中）
    POSITIVE_SELECTION_POSITIONS = [2, 3, 4, 6, 7, 8, 10, 12, 13, 15, 17, 20, 
                                    21, 22, 23, 25, 26, 28, 29, 30, 31, 32, 34, 
                                    35, 36, 44, 46, 49, 57, 75, 96]
    
    # 分析参数
    N_AA = 100  # 分析N端100个氨基酸
    
    # 创建输出目录
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    
    # 初始化分析器
    analyzer = EnhancedPositionImportanceAnalyzer(
        model_name="facebook/esm2_t33_650M_UR50D",
        device="cpu"  # 或 "cuda"
    )
    
    try:
        # 1. 加载数据（简化版）
        print("\n" + "=" * 60)
        print("STEP 1: Loading Data")
        print("=" * 60)
        
        import pandas as pd
        
        # 加载聚类文件
        cluster_df = pd.read_csv(CLUSTER_FILE)
        gene_cluster_dict = {}
        for _, row in cluster_df.iterrows():
            genes = str(row['Gene_List']).split(';')
            for gene in genes:
                gene_cluster_dict[gene.strip()] = str(row['Cluster'])
        
        # 加载FASTA文件
        sequences_dict = {}
        with open(FASTA_FILE, 'r') as f:
            current_gene = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_gene and current_seq:
                        sequences_dict[current_gene] = ''.join(current_seq)
                    
                    header = line[1:]
                    if '|' in header:
                        current_gene = header.split('|')[0]
                    else:
                        current_gene = header.split()[0]
                    
                    current_seq = []
                elif line:
                    current_seq.append(line)
            
            if current_gene and current_seq:
                sequences_dict[current_gene] = ''.join(current_seq)
        
        # 对齐数据
        common_genes = set(sequences_dict.keys()) & set(gene_cluster_dict.keys())
        sequences_dict = {g: sequences_dict[g] for g in common_genes}
        gene_cluster_dict = {g: gene_cluster_dict[g] for g in common_genes}
        
        print(f"Loaded {len(common_genes)} genes with both sequence and cluster labels")
        
        # 2. 准备N端序列数据
        print("\n" + "=" * 60)
        print("STEP 2: Preparing N-terminal Sequences")
        print("=" * 60)
        
        n_term_sequences = []
        labels = []
        
        for gene, seq in sequences_dict.items():
            if len(seq) >= N_AA:
                n_term_sequences.append(seq[:N_AA])
                labels.append(gene_cluster_dict[gene])
        
        print(f"Prepared {len(n_term_sequences)} N-terminal {N_AA}AA sequences")
        print(f"Number of classes: {len(set(labels))}")
        
        # 打印类别分布
        class_counts = pd.Series(labels).value_counts()
        print(f"Class distribution: {class_counts.to_dict()}")
        print(f"Top 5 classes: {class_counts.head(5).to_dict()}")
        
        # 3. 运行6次重复分析
        print("\n" + "=" * 60)
        print("STEP 3: Running 6-Repeat Position Importance Analysis")
        print("=" * 60)
        
        all_results = analyzer.six_repeats_analysis(n_term_sequences, labels, N_AA)
        
        # 4. 与正选择位点比较
        print("\n" + "=" * 60)
        print("STEP 4: Comparing with Positive Selection Positions")
        print("=" * 60)
        
        key_positions = all_results['integrated']['position_importance'].get('key_positions', [])
        comparison_results = analyzer.compare_with_positive_selection(
            key_positions, POSITIVE_SELECTION_POSITIONS, N_AA
        )
        
        # 5. 可视化
        print("\n" + "=" * 60)
        print("STEP 5: Generating Visualizations")
        print("=" * 60)
        
        analyzer.visualize_integrated_results(
            all_results['integrated'], 
            POSITIVE_SELECTION_POSITIONS,
            OUTPUT_DIR
        )
        
        # 6. 生成综合报告
        print("\n" + "=" * 60)
        print("STEP 6: Generating Comprehensive Report")
        print("=" * 60)
        
        analyzer.generate_comprehensive_report(
            all_results['integrated'],
            comparison_results,
            f"{OUTPUT_DIR}/comprehensive_analysis_report.txt"
        )
        
        # 7. 保存详细数据（修复JSON序列化问题）
        print("\n" + "=" * 60)
        print("STEP 7: Saving Detailed Data (Fixed JSON Serialization)")
        print("=" * 60)
        
        # 保存整合的位置重要性数据
        importance_df = pd.DataFrame({
            'position': range(1, N_AA + 1),
            'mean_importance': all_results['integrated']['position_importance']['mean_across_repeats'],
            'std_importance': all_results['integrated']['position_importance']['std_across_repeats'],
            'median_importance': all_results['integrated']['position_importance']['median_across_repeats'],
            'coefficient_of_variation': all_results['integrated']['position_importance']['coefficient_of_variation'],
            'is_key_position': [1 if (i+1) in key_positions else 0 for i in range(N_AA)],
            'is_positive_selection': [1 if (i+1) in POSITIVE_SELECTION_POSITIONS else 0 for i in range(N_AA)]
        })
        importance_df.to_csv(f"{OUTPUT_DIR}/integrated_position_importance.csv", index=False)
        
        # 保存比较结果（使用可序列化转换）
        comparison_serializable = convert_to_serializable(comparison_results)
        with open(f"{OUTPUT_DIR}/comparison_results.json", 'w') as f:
            json.dump(comparison_serializable, f, indent=2, cls=NumpyEncoder)
        
        # 保存元数据（使用可序列化转换）
        metadata = {
            'analysis_date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
            'n_sequences': len(n_term_sequences),
            'n_classes': len(set(labels)),
            'n_aa_analyzed': N_AA,
            'positive_selection_positions': POSITIVE_SELECTION_POSITIONS,
            'key_positions_identified': key_positions,
            'accuracy_summary': all_results['integrated']['accuracy'],
            'stability_summary': all_results['integrated']['stability'],
            'n_repeats_successful': all_results['n_successful_repeats']
        }
        
        metadata_serializable = convert_to_serializable(metadata)
        with open(f"{OUTPUT_DIR}/analysis_metadata.json", 'w') as f:
            json.dump(metadata_serializable, f, indent=2, cls=NumpyEncoder)
        
        print(f"\n✅ ANALYSIS COMPLETED SUCCESSFULLY!")
        print(f"📁 Results saved to: {OUTPUT_DIR}/")
        print(f"\n📊 KEY OUTPUT FILES:")
        print(f"   1. {OUTPUT_DIR}/comprehensive_analysis_report.txt - Detailed analysis report")
        print(f"   2. {OUTPUT_DIR}/integrated_position_importance.csv - Position importance scores")
        print(f"   3. {OUTPUT_DIR}/comparison_results.json - SVM vs positive selection comparison")
        print(f"   4. {OUTPUT_DIR}/plots/ - All visualization charts")
        print(f"   5. {OUTPUT_DIR}/analysis_metadata.json - Analysis metadata")
        print(f"   6. repeat_*/ - Individual repeat results (6 folders)")
        
        # 打印关键发现
        print(f"\n🔬 KEY FINDINGS:")
        print(f"   - Mean accuracy across 6 repeats: {all_results['integrated']['accuracy']['mean_accuracy']:.4f}")
        print(f"   - Number of key positions identified: {len(key_positions)}")
        print(f"   - Overlap with positive selection positions: {comparison_results['n_overlap']}")
        print(f"   - Statistical significance of overlap: p = {comparison_results['p_value']:.6f}")
        
        if comparison_results['p_value'] < 0.05:
            print(f"   ✅ Significant overlap between SVM key positions and positive selection sites!")
            print(f"   Overlap positions: {sorted(comparison_results['overlap_positions'])}")
        else:
            print(f"   ⚠️ No significant overlap found")
        
        print(f"\n🧪 RECOMMENDATIONS:")
        print(f"   1. Validate top key positions experimentally")
        if comparison_results['p_value'] < 0.05 and comparison_results['n_overlap'] > 0:
            print(f"   2. Focus on positions with both SVM importance and positive selection: {sorted(comparison_results['overlap_positions'])}")
        print(f"   3. Consider structural context of identified positions")
        
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()