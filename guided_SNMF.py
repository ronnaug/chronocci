"""
Guided Spatial Non-negative Matrix Factorization (Guided-SNMF).
Tailored for ALCL with adaptive expression gating and scikit-learn API.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
import logging

# Настройка логирования для биоинформатического пайплайна
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class GuidedSNMF(nn.Module):
    def __init__(self, n_genes, n_spots, n_factors, cd30_idx, jak_stat_indices, t_cell_indices):
        super(GuidedSNMF, self).__init__()
        self.n_factors = n_factors
        
        # Случайный старт в лог-пространстве (для стабильности softplus)
        W_init = torch.randn(n_genes, n_factors) * 0.02 - 1.0 
        H_init = torch.randn(n_factors, n_spots) * 0.02 - 1.0
        
        # Закладываем эмбриональный профиль опухоли в Фактор 0
        W_init[cd30_idx, 0] = 2.0          
        W_init[jak_stat_indices, 0] = 1.5  
        W_init[t_cell_indices, 0] = -15.0  # Глубокий минус в softplus даст чистый 0
        
        self.W_raw = nn.Parameter(W_init)
        self.H_raw = nn.Parameter(H_init)
        
    @property
    def W(self):
        return F.softplus(self.W_raw)

    @property
    def H(self):
        return F.softplus(self.H_raw)


class UniversalGuidedSpatialNMFLoss(nn.Module):
    def __init__(self, cd30_idx, cd45_idx, jak_stat_indices, t_cell_indices,
                 lambda_spatial=1.5, lambda_cd30=30.0, lambda_jak=20.0, 
                 lambda_t=40.0, lambda_background=100.0, eps=1e-8):
        super(UniversalGuidedSpatialNMFLoss, self).__init__()
        self.cd30_idx = cd30_idx
        self.cd45_idx = cd45_idx
        self.jak_stat_indices = jak_stat_indices
        self.t_cell_indices = t_cell_indices
        
        self.l_spatial = lambda_spatial
        self.l_cd30 = lambda_cd30
        self.l_jak = lambda_jak
        self.l_t = lambda_t
        self.l_background = lambda_background
        self.eps = eps

    def forward(self, X_true, W, H, S_matrix):
        X_pred = torch.mm(W, H) 
        loss_recon = F.mse_loss(X_pred, X_true)
        
        H_neighbors = torch.mm(H, S_matrix.t())
        loss_spatial = F.mse_loss(H, H_neighbors)
        
        W_norm = W / (W.sum(dim=0, keepdims=True) + self.eps)
        
        # Динамический поиск опухолевого фактора
        tumor_scores = W_norm[self.cd30_idx, :] + W_norm[self.jak_stat_indices, :].mean(dim=0)
        tumor_k = torch.argmax(tumor_scores).item() 
        
        # Двухэтапный относительный гейтинг
        total_hvg_expression_per_spot = X_true.sum(dim=0) + self.eps
        cd30_share_per_spot = X_true[self.cd30_idx, :] / total_hvg_expression_per_spot
        top_cd30_shares = torch.topk(cd30_share_per_spot, k=min(10, cd30_share_per_spot.shape[0])).values
        mean_top_cd30_share = top_cd30_shares.mean()
        
        gate_weight = torch.clamp((mean_top_cd30_share - 0.01) / 0.01, min=0.0, max=1.0)
        
        loss_cd30 = gate_weight * (torch.relu(0.05 - W_norm[self.cd30_idx, tumor_k])**2)
        loss_jak = gate_weight * (torch.relu(0.02 - torch.mean(W_norm[jak_stat_indices, tumor_k]))**2)
        loss_t_identity = torch.sum(W_norm[self.t_cell_indices, tumor_k])**2
        
        if self.cd45_idx is not None:
            loss_cd45 = torch.relu(0.01 - W_norm[self.cd45_idx, tumor_k])**2
        else:
            loss_cd45 = 0.0
            
        loss_zero_tissue = (1.0 - gate_weight) * torch.mean(H[tumor_k, :])**2
        
        total_loss = (loss_recon + 
                      self.l_spatial * loss_spatial + 
                      self.l_cd30 * loss_cd30 + 
                      self.l_jak * loss_jak + 
                      self.l_t * (loss_t_identity + loss_cd45) +
                      self.l_background * loss_zero_tissue)
                      
        return total_loss, loss_recon.item(), tumor_k


# =========================================================================
# 🚀 ВЫСОКОУРОВНЕВЫЙ ИНТЕРФЕЙС ПАКЕТА (Scikit-Learn Style API)
# =========================================================================
class SpatialGuidedDeconvolution:
    """
    Основной интерфейс деконволюции. Автоматически парсит маркеры.
    По умолчанию настроен на биологические сигнатуры АККЛ (ALCL).
    """
    def __init__(self, 
                 n_factors=15, 
                 main_onco_marker="TNFRSF8", # CD30
                 pan_leukocyte_marker="PTPRC", # CD45
                 pathway_signature=None, 
                 t_cell_canonical=None,
                 lambda_spatial=1.5,
                 lambda_onco=30.0,
                 lambda_pathway=20.0,
                 lambda_t_loss=40.0,
                 lambda_bg_kill=100.0):
        
        self.n_factors = n_factors
        self.main_marker = main_onco_marker
        self.pan_marker = pan_leukocyte_marker
        
        # Зашиваем дефолтные маркеры АККЛ, если пользователь не передал свои
        self.pathway_sig = pathway_signature if pathway_signature is not None else ["STAT3", "SOCS3", "PIM1"]
        self.t_canon = t_cell_canonical if t_cell_canonical is not None else ["CD3D", "CD3E", "CD5", "CD28"]
        
        # Коэффициенты регуляризации
        self.l_spatial = lambda_spatial
        self.l_onco = lambda_onco
        self.l_pathway = lambda_pathway
        self.l_t = lambda_t_loss
        self.l_bg = lambda_bg_kill
        
        # Внутренние переменные для результатов
        self.W_ = None
        self.H_ = None
        self.tumor_factor_index_ = None

    def _parse_gene_indices(self, var_names):
        """Внутренний хелпер для безопасного перевода названий генов в индексы"""
        gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
        
        if self.main_marker not in gene_to_idx:
            raise ValueError(f"Критический онкомаркер '{self.main_marker}' отсутствует в var_names датасета!")
            
        cd30_idx = gene_to_idx[self.main_marker]
        cd45_idx = gene_to_idx.get(self.pan_marker, None)
        
        # Парсим списки генов, отсекая те, которых нет в текущей матрице HVG
        jak_indices = [gene_to_idx[g] for g in self.pathway_sig if g in gene_to_idx]
        t_indices = [gene_to_idx[g] for g in self.t_canon if g in gene_to_idx]
        
        if not jak_indices:
            logging.warning("Ни один ген из сигнатуры пути не найден в var_names. Регуляризация путей будет ослаблена.")
        if not t_indices:
            logging.warning("Канонические Т-маркеры отсутствуют в var_names. Контроль Т-идентичности отключен.")
            
        return cd30_idx, cd45_idx, jak_indices, t_indices

    def fit(self, adata_hvg, S_matrix, n_epochs=1000, lr=0.01, device='cpu'):
        """
        Обучает модель деконволюции на основе переданного AnnData объекта и матрицы смежности.
        """
        logging.info("Парсинг биологических маркеров и подготовка тензоров...")
        cd30_idx, cd45_idx, jak_indices, t_indices = self._parse_gene_indices(adata_hvg.var_names)
        
        # Конвертируем данные в тензоры [Genes x Spots]
        X_raw = adata_hvg.X.toarray().T if hasattr(adata_hvg.X, "toarray") else adata_hvg.X.copy().T
        X_tensor = torch.tensor(X_raw, dtype=torch.float32).to(device)
        S_tensor = torch.tensor(S_matrix, dtype=torch.float32).to(device)
        
        n_genes, n_spots = X_tensor.shape
        
        # Инициализируем модель
        model = GuidedSNMF(
            n_genes=n_genes, n_spots=n_spots, n_factors=self.n_factors,
            cd30_idx=cd30_idx, jak_stat_indices=jak_indices, t_cell_indices=t_indices
        ).to(device)
        
        criterion = UniversalGuidedSpatialNMFLoss(
            cd30_idx=cd30_idx, cd45_idx=cd45_idx, jak_stat_indices=jak_indices, t_cell_indices=t_indices,
            lambda_spatial=self.l_spatial, lambda_cd30=self.l_onco, lambda_jak=self.l_pathway,
            lambda_t=self.l_t, lambda_background=self.l_bg
        )
        
        optimizer = optim.Adam(model.parameters(), lr=lr)
        
        logging.info(f"Старт оптимизации Guided Spatial NMF на {device} ({n_epochs} эпох)...")
        for epoch in range(1, n_epochs + 1):
            optimizer.zero_grad()
            total_loss, _, active_tumor_idx = criterion(X_tensor, model.W, model.H, S_tensor)
            total_loss.backward()
            optimizer.step()
            
            if epoch % 200 == 0 or epoch == 1:
                logging.info(f"  Эпоха {epoch:04d}/{n_epochs} | Loss: {total_loss.item():.4f} | Опухолевый фактор: {active_tumor_idx}")
                
        # Фиксируем результаты обучения внутри инстанса класса
        self.W_ = model.W.detach().cpu().numpy()
        self.H_ = model.H.detach().cpu().numpy()
        self.tumor_factor_index_ = active_tumor_idx
        logging.info("Обучение успешно завершено.")
        return self

    def get_tumor_weights(self):
        """Возвращает изолированный одномерный вектор пространственной активности опухоли"""
        if self.H_ is None:
            raise ValueError("Модель еще не обучена! Сначала вызовите метод .fit()")
        return self.H_[self.tumor_factor_index_, :]