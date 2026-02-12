"""Custom scoring functions for voxel connectivity models.

This module provides hybrid scoring that combines voxel-level and region-level
mean squared relative error metrics for evaluating connectivity model predictions.
"""
from __future__ import division
import numpy as np

from sklearn.metrics import make_scorer
from sklearn.metrics._regression import _check_reg_targets

from mcmodels.core import Mask
from mcmodels.utils import squared_norm


class HybridScorer(object):
    """Hybrid scorer combining voxel and regional scoring metrics.
    
    Attributes:
        DEFAULT_STRUCTURE_SET_ID: Default Allen CCF structure set identifier
        cache: VoxelModelCache instance
        structure_ids: List of structure IDs to use for scoring
    """

    DEFAULT_STRUCTURE_SET_ID = 687527945

    @staticmethod
    def voxel_scorer():
        """Create a voxel-level scorer.
        
        Returns:
            sklearn scorer: Mean squared relative error scorer for voxel predictions
        """
        return make_scorer(mean_squared_relative_error, greater_is_better=False)

    @staticmethod
    def regional_scorer(**kwargs):
        """Create a region-level scorer.
        
        Args:
            **kwargs: Additional arguments for regional scoring
            
        Returns:
            sklearn scorer: Regional mean squared relative error scorer
        """
        return make_scorer(regional_mean_squared_relative_error, greater_is_better=False, **kwargs)

    @property
    def _default_structure_ids(self):
        """Get default structure IDs from Allen CCF.
        
        Returns:
            list: Structure IDs excluding specific structures (934, 1009)
        """
        structure_tree = self.cache.get_structure_tree()
        structures = structure_tree.get_structures_by_set_id([self.DEFAULT_STRUCTURE_SET_ID])

        return [s['id'] for s in structures if s['id'] not in (934, 1009)]

    def __init__(self, cache, structure_ids=None):
        """Initialize HybridScorer.
        
        Args:
            cache: VoxelModelCache instance
            structure_ids: Optional list of structure IDs (uses default if None)
        """
        self.cache = cache
        self.structure_ids = structure_ids

        if self.structure_ids is None:
            self.structure_ids = self._default_structure_ids

    @property
    def scoring_dict(self):
        """Create dictionary of scoring functions.
        
        Returns:
            dict: Dictionary with 'voxel' and 'regional' scorers
        """
        def get_nnz_assigned(key):
            """Get non-zero assigned region indices from key.
            
            Args:
                key: Region assignment key array
                
            Returns:
                np.ndarray: Non-zero region indices
            """
            assigned = np.unique(key)
            if assigned[0] == 0:
                return assigned[1:]
            return assigned

        # Target is whole brain
        target_mask = Mask.from_cache(cache=self.cache, hemisphere_id=3)

        ipsi_key = target_mask.get_key(structure_ids=self.structure_ids, hemisphere_id=2)
        contra_key = target_mask.get_key(structure_ids=self.structure_ids, hemisphere_id=1)

        reg_kwargs = dict(ipsi_key=ipsi_key,
                          contra_key=contra_key,
                          ipsi_regions=get_nnz_assigned(ipsi_key),
                          contra_regions=get_nnz_assigned(contra_key))

        return dict(voxel=self.voxel_scorer(), regional=self.regional_scorer(**reg_kwargs))


def log_mean_squared_relative_error(y_true, y_pred):
    """Calculate log-transformed mean squared relative error.
    
    Args:
        y_true: True values
        y_pred: Predicted values
        
    Returns:
        float: Log-transformed MSRE
    """
    log = lambda x: np.log10(x + 1e-8)
    return mean_squared_relative_error(log(y_true), log(y_pred))


def log_regional_mean_squared_relative_error(y_true, y_pred, **kwargs):
    """Calculate log-transformed regional mean squared relative error.
    
    Args:
        y_true: True values
        y_pred: Predicted values
        **kwargs: Additional arguments for regional scoring
        
    Returns:
        float: Log-transformed regional MSRE
    """
    log = lambda x: np.log10(x + 1e-8)
    return regional_mean_squared_relative_error(log(y_true), log(y_pred), **kwargs)


class LogHybridScorer(HybridScorer):
    """Hybrid scorer using log-transformed metrics.
    
    Inherits from HybridScorer but uses log-transformed scoring functions.
    """

    @staticmethod
    def voxel_scorer():
        """Create a log-transformed voxel-level scorer.
        
        Returns:
            sklearn scorer: Log-transformed MSRE scorer for voxel predictions
        """
        return make_scorer(log_mean_squared_relative_error, greater_is_better=False)

    @staticmethod
    def regional_scorer(**kwargs):
        """Create a log-transformed region-level scorer.
        
        Args:
            **kwargs: Additional arguments for regional scoring
            
        Returns:
            sklearn scorer: Log-transformed regional MSRE scorer
        """
        return make_scorer(log_regional_mean_squared_relative_error,
                           greater_is_better=False, **kwargs)


def unionize(v, ipsi_key, contra_key, ipsi_regions, contra_regions):
    """Unionize voxel-level values to regional values.
    
    Aggregates voxel predictions to regional predictions by summing values
    for each region defined in ipsi_key and contra_key.
    
    Args:
        v: Voxel values array of shape (n_samples, n_voxels)
        ipsi_key: Ipsilateral region assignment key
        contra_key: Contralateral region assignment key
        ipsi_regions: Ipsilateral region indices
        contra_regions: Contralateral region indices
        
    Returns:
        np.ndarray: Regional values of shape (n_samples, n_regions)
        
    Raises:
        ValueError: If keys are incompatible or v has wrong number of columns
    """
    if ipsi_key.shape != contra_key.shape:
        raise ValueError("Keys are incompatible: ipsi and contra keys must have same shape")

    v = np.atleast_2d(v)
    if v.shape[1] != ipsi_key.size:  # or contra, doesn't matter
        raise ValueError("Key must be the same size as the n columns in vector!")

    j = 0
    result = np.empty((v.shape[0], len(ipsi_regions) + len(contra_regions)))
    for key, regions in zip((ipsi_key, contra_key), (ipsi_regions, contra_regions)):
        for k in regions:
            result[:, j] = v[:, np.where(key == k)[0]].sum(axis=1)
            j += 1

    return result


def mean_squared_relative_error(y_true, y_pred, multioutput='uniform_average'):
    """Calculate mean squared relative error.
    
    Computes: 2 * ||y_true - y_pred||^2 / (||y_true||^2 + ||y_pred||^2)
    
    Args:
        y_true: True values
        y_pred: Predicted values
        multioutput: How to aggregate multioutput results
        
    Returns:
        float: Mean squared relative error
    """
    _, y_true, y_pred, _ = _check_reg_targets(y_true, y_pred, multioutput)
    result = 2 * squared_norm(y_true - y_pred) / (squared_norm(y_true) + squared_norm(y_pred))
    return result


def regional_mean_squared_relative_error(y_true, y_pred, **kwargs):
    """Calculate regional mean squared relative error.
    
    First unionizes voxel predictions to regional predictions, then computes MSRE.
    
    Args:
        y_true: True voxel values
        y_pred: Predicted voxel values
        **kwargs: Must include ipsi_key, contra_key, ipsi_regions, contra_regions
        
    Returns:
        float: Regional mean squared relative error
        
    Raises:
        ValueError: If required kwargs are missing
    """
    try:
        ipsi_key = kwargs.pop('ipsi_key')
        contra_key = kwargs.pop('contra_key')
        ipsi_regions = kwargs.pop('ipsi_regions')
        contra_regions = kwargs.pop('contra_regions')
    except KeyError:
        raise ValueError("Must be called with 'ipsi_key', 'contra_key', 'ipsi_regions', "
                         "and 'contra_regions' kwargs")

    y_true = unionize(y_true, ipsi_key, contra_key, ipsi_regions, contra_regions)
    y_pred = unionize(y_pred, ipsi_key, contra_key, ipsi_regions, contra_regions)

    return mean_squared_relative_error(y_true, y_pred, **kwargs)


def mse_rel():
    """Create a mean squared relative error scorer.
    
    Returns:
        sklearn scorer: MSRE scorer
    """
    return make_scorer(mean_squared_relative_error, greater_is_better=False)


def regional_mse_rel(**kwargs):
    """Create a regional mean squared relative error scorer.
    
    Args:
        **kwargs: Additional arguments for regional scoring
        
    Returns:
        sklearn scorer: Regional MSRE scorer
    """
    return make_scorer(regional_mean_squared_relative_error, greater_is_better=False, **kwargs)
