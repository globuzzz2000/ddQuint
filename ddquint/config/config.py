#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration module for the ddQuint pipeline with standard deviation-based tolerances.

This module provides comprehensive configuration management for:
1. Clustering parameters and algorithm settings
2. Expected centroid definitions for up to 10 chromosomes
3. Standard deviation-based copy number classification with tolerance multiplier
4. Visualization settings and color schemes
5. File management and template parsing options

The Config class implements a singleton pattern to ensure consistent
settings across all pipeline modules.
"""

import os
import json
import logging
from multiprocessing import cpu_count
from typing import Dict, List, Any

from ..config.exceptions import ConfigError

logger = logging.getLogger(__name__)

class Config:
    """
    Central configuration settings for the ddQuint pipeline with singleton pattern.
    
    This class provides configuration management for clustering parameters,
    visualization settings, copy number thresholds, and file management.
    Implements singleton pattern to ensure consistent settings across modules.
    
    Attributes:
        DEBUG_MODE: Enable debug logging mode
        EXPECTED_CENTROIDS: Target centroids for clustering
        HDBSCAN_MIN_CLUSTER_SIZE: Minimum cluster size for HDBSCAN
        
    Example:
        >>> config = Config.get_instance()
        >>> config.DEBUG_MODE = True
        >>> chroms = config.get_chromosome_keys()
    """
    
    # Singleton instance
    _instance = None

    #############################################################################
    #                           Expected Centroids
    #############################################################################
    # Define expected centroids for targets (maximum 10 chromosomes)
    # Format: { "target_name": [Ch1Amplitude, Ch2Amplitude] }
    EXPECTED_CENTROIDS = {
        "Negative": [1000, 900],
        "Chrom1":   [1000, 2300],
        "Chrom2":   [1800, 2200],
        "Chrom3":   [2400, 1750],
        "Chrom4":   [3100, 1300],
        "Chrom5":   [3500, 900]
    }
    
    # Tolerance for matching clusters to targets
    BASE_TARGET_TOLERANCE = 750
    
    # Scale factor limits for adaptive tolerance
    SCALE_FACTOR_MIN = 0.5
    SCALE_FACTOR_MAX = 1.0
    

    #############################################################################
    #                           Clustering Settings
    #############################################################################
    # HDBSCAN clustering parameters
    HDBSCAN_MIN_CLUSTER_SIZE = 4
    HDBSCAN_MIN_SAMPLES = 70
    HDBSCAN_EPSILON = 0.06
    HDBSCAN_METRIC = 'euclidean'
    HDBSCAN_CLUSTER_SELECTION_METHOD = 'eom'
    
    # Minimum data points required for clustering
    MIN_POINTS_FOR_CLUSTERING = 50
    

    #############################################################################
    #                    Standard Deviation-Based Copy Number Settings
    #############################################################################
    # Copy number calculation parameters
    MIN_USABLE_DROPLETS = 3000
    COPY_NUMBER_MEDIAN_DEVIATION_THRESHOLD = 0.15  # 15% deviation threshold
    COPY_NUMBER_BASELINE_MIN_CHROMS = 3  # Minimum chromosomes for baseline calc
    
    # Expected copy number values for each chromosome (baseline for calculations)
    EXPECTED_COPY_NUMBERS = {
        "Chrom1": 0.9716,
        "Chrom2": 1.0052,
        "Chrom3": 1.0278,
        "Chrom4": 0.9912,
        "Chrom5": 1.0035
    }
    
    # Standard deviation for each chromosome (empirically determined)
    EXPECTED_STANDARD_DEVIATION = {
        "Chrom1": 0.0312,
        "Chrom2": 0.0241,
        "Chrom3": 0.0290,
        "Chrom4": 0.0242,
        "Chrom5": 0.0230
    }
    
    # Tolerance multiplier for standard deviation-based classification
    TOLERANCE_MULTIPLIER = 3
    
    # Aneuploidy target multipliers (multiplicative factors for expected values)
    ANEUPLOIDY_TARGETS = {
        "low": 0.75,   # Deletion target (expected * 0.75)
        "high": 1.25   # Duplication target (expected * 1.25)
    }

    #############################################################################
    #                           Visualization Settings
    #############################################################################
    # Plot dimensions and settings
    COMPOSITE_FIGURE_SIZE = (16, 11)
    INDIVIDUAL_FIGURE_SIZE = (6, 5)
    COMPOSITE_PLOT_SIZE = (5, 5)
    
    # DPI settings for different plot types
    INDIVIDUAL_PLOT_DPI = 300      # High resolution for standalone plots
    COMPOSITE_PLOT_DPI = 200       # Medium resolution for composite images
    PLACEHOLDER_PLOT_DPI = 150     # Lower resolution for placeholder plots
    
    # Axis limits
    X_AXIS_MIN = 0
    X_AXIS_MAX = 3000
    Y_AXIS_MIN = 0
    Y_AXIS_MAX = 5000
    
    # Grid settings
    X_GRID_INTERVAL = 500
    Y_GRID_INTERVAL = 1000
    
    # Color scheme for targets (up to 10 chromosomes)
    DEFAULT_COLOR_PALETTE = [
        "#1f77b4", "#f59a23", "#7ec638", "#16d9ff", "#f65352",
        "#82218b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    ]
    SPECIAL_COLOR_DEFAULTS = {
        "Negative": "#1f77b4",
        "Unknown":  "#c7c7c7",
    }

    # Copy number state highlighting colors
    ANEUPLOIDY_FILL_COLOR = "#E6B8E6"  # Light purple (for definitive aneuploidies)
    ANEUPLOIDY_VALUE_FILL_COLOR = "#D070D0"  # Darker purple (for aneuploidy values)
    BUFFER_ZONE_FILL_COLOR = "#B0B0B0"  # Dark grey (for buffer zone samples - entire row)
    BUFFER_ZONE_VALUE_FILL_COLOR = "#808080"  # Darker grey (for buffer zone values - not used now)
    
    #############################################################################
    #                           File Management
    #############################################################################
    # Default output directories
    GRAPHS_DIR_NAME = "Graphs"
    COMPOSITE_IMAGE_FILENAME = "Graph_Overview.png"
    CSV_EXTENSION = '.csv'        # File name patterns
    
    #############################################################################
    #                           Template Parsing
    #############################################################################
    # Template search parameters
    TEMPLATE_SEARCH_PARENT_LEVELS = 2  # How many parent directories to search up
    TEMPLATE_PATTERN = "{dir_name}.csv"  # Template file naming pattern
    
    #############################################################################
    #                           Well Management
    #############################################################################
    # 96-well plate layout
    PLATE_ROWS = list('ABCDEFGH')
    PLATE_COLS = [str(i) for i in range(1, 13)]
    WELL_FORMAT = "{row}{col:02d}"  # e.g., "A01"
    
    def __init__(self):
        """Initialize Config instance with default values."""
        pass
    
    @classmethod
    def get_instance(cls) -> 'Config':
        """
        Get the singleton instance of Config.
        
        Returns:
            Config: Singleton instance
        """
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance
    
    @classmethod
    def get_chromosome_keys(cls) -> List[str]:
        """
        Get all chromosome keys from expected centroids.
        
        Returns:
            List of chromosome keys sorted numerically
            
        Example:
            >>> config = Config.get_instance()
            >>> chroms = config.get_chromosome_keys()
            >>> chroms
            ['Chrom1', 'Chrom2', 'Chrom3']
        """
        return sorted([key for key in cls.EXPECTED_CENTROIDS.keys() 
                      if key.startswith('Chrom')], 
                     key=lambda x: int(x.replace('Chrom', '')))
    
    @classmethod
    def get_ordered_labels(cls) -> List[str]:
        """
        Get ordered labels including all chromosomes.
        
        Returns:
            List of labels in processing order
        """
        return ['Negative'] + cls.get_chromosome_keys() + ['Unknown']
    
    @classmethod
    def get_tolerance_for_chromosome(cls, chrom_name: str) -> float:
        """
        Get the tolerance value for a specific chromosome based on its standard deviation.
        
        Args:
            chrom_name: Chromosome name (e.g., 'Chrom1')
            
        Returns:
            Tolerance value calculated as std_dev * multiplier
            
        Raises:
            ConfigError: If chromosome not found in configuration
        """
        if chrom_name not in cls.EXPECTED_STANDARD_DEVIATION:
            error_msg = f"Unknown chromosome for standard deviation: {chrom_name}"
            logger.error(error_msg)
            raise ConfigError(error_msg, config_key="EXPECTED_STANDARD_DEVIATION")
        
        std_dev = cls.EXPECTED_STANDARD_DEVIATION[chrom_name]
        tolerance = std_dev * cls.TOLERANCE_MULTIPLIER
        
        logger.debug(f"{chrom_name}: std_dev={std_dev:.4f}, tolerance={tolerance:.4f}")
        return tolerance
    
    @classmethod
    def classify_copy_number_state(cls, chrom_name: str, copy_number: float) -> str:
        """
        Classify a copy number value using standard deviation-based tolerances.
        
        Uses chromosome-specific standard deviations with a tolerance multiplier to
        determine the classification state for copy number analysis.
        
        Args:
            chrom_name: Chromosome name (e.g., 'Chrom1')
            copy_number: Copy number value to classify
            
        Returns:
            Classification string: 'euploid', 'buffer_zone', or 'aneuploidy'
            
        Raises:
            ConfigError: If chromosome not found in configuration
            
        Example:
            >>> config = Config.get_instance()
            >>> state = config.classify_copy_number_state('Chrom1', 1.0)
            >>> state
            'euploid'
        """
        if chrom_name not in cls.EXPECTED_COPY_NUMBERS:
            error_msg = f"Unknown chromosome: {chrom_name}"
            logger.error(error_msg)
            raise ConfigError(error_msg, config_key="EXPECTED_COPY_NUMBERS")
        
        # Get expected value and tolerance for this chromosome
        expected = cls.EXPECTED_COPY_NUMBERS[chrom_name]
        tolerance = cls.get_tolerance_for_chromosome(chrom_name)
        
        # Define euploid range using chromosome-specific tolerance
        euploid_min = expected - tolerance
        euploid_max = expected + tolerance
        
        # Define aneuploidy target ranges using the same tolerance
        deletion_target = expected * cls.ANEUPLOIDY_TARGETS["low"]
        duplication_target = expected * cls.ANEUPLOIDY_TARGETS["high"]
        
        deletion_min = deletion_target - tolerance
        deletion_max = deletion_target + tolerance
        duplication_min = duplication_target - tolerance
        duplication_max = duplication_target + tolerance
        
        logger.debug(f"{chrom_name} classification ranges:")
        logger.debug(f"  Euploid: [{euploid_min:.4f}, {euploid_max:.4f}]")
        logger.debug(f"  Deletion: [{deletion_min:.4f}, {deletion_max:.4f}]")
        logger.debug(f"  Duplication: [{duplication_min:.4f}, {duplication_max:.4f}]")
        logger.debug(f"  Copy number: {copy_number:.4f}")
        
        # Check if in euploid range
        if euploid_min <= copy_number <= euploid_max:
            logger.debug(f"  -> euploid")
            return 'euploid'
        
        # Check if in aneuploidy ranges
        if (deletion_min <= copy_number <= deletion_max or 
            duplication_min <= copy_number <= duplication_max):
            logger.debug(f"  -> aneuploidy")
            return 'aneuploidy'
        
        # Otherwise, it's in the buffer zone
        logger.debug(f"  -> buffer_zone")
        return 'buffer_zone'
    
    @classmethod
    def get_copy_number_ranges(cls, chrom_name: str) -> Dict[str, tuple]:
        """
        Get copy number ranges for a specific chromosome using standard deviation-based tolerances.
        
        Args:
            chrom_name: Chromosome name (e.g., 'Chrom1')
            
        Returns:
            Dictionary with ranges for each classification
            
        Raises:
            ConfigError: If chromosome not found in configuration
        """
        if chrom_name not in cls.EXPECTED_COPY_NUMBERS:
            error_msg = f"Unknown chromosome: {chrom_name}"
            logger.error(error_msg)
            raise ConfigError(error_msg, config_key="EXPECTED_COPY_NUMBERS")
        
        expected = cls.EXPECTED_COPY_NUMBERS[chrom_name]
        tolerance = cls.get_tolerance_for_chromosome(chrom_name)
        
        # Calculate ranges using chromosome-specific tolerance
        euploid_range = (
            expected - tolerance,
            expected + tolerance
        )
        
        # Use multiplicative targets for aneuploidy
        deletion_target = expected * cls.ANEUPLOIDY_TARGETS["low"]
        duplication_target = expected * cls.ANEUPLOIDY_TARGETS["high"]
        
        deletion_range = (
            deletion_target - tolerance,
            deletion_target + tolerance
        )
        
        duplication_range = (
            duplication_target - tolerance,
            duplication_target + tolerance
        )
        
        return {
            'euploid': euploid_range,
            'deletion': deletion_range,
            'duplication': duplication_range
        }
    
    @classmethod
    def get_plot_dpi(cls, plot_type: str = 'individual') -> int:
        """
        Get DPI setting for different plot types.
        
        Args:
            plot_type: Type of plot ('individual', 'composite', or 'placeholder')
            
        Returns:
            DPI value for the specified plot type
            
        Raises:
            ConfigError: If unknown plot type is specified
            
        Example:
            >>> config = Config.get_instance()
            >>> dpi = config.get_plot_dpi('composite')
            >>> dpi
            200
        """
        dpi_mapping = {
            'individual': cls.INDIVIDUAL_PLOT_DPI,
            'composite': cls.COMPOSITE_PLOT_DPI,
            'placeholder': cls.PLACEHOLDER_PLOT_DPI
        }
        
        if plot_type not in dpi_mapping:
            error_msg = f"Unknown plot type: {plot_type}. Valid types: {list(dpi_mapping.keys())}"
            logger.error(error_msg)
            raise ConfigError(error_msg, config_key="plot_type")
        
        return dpi_mapping[plot_type]
    
    @classmethod
    def load_from_file(cls, filepath: str) -> bool:
        """
        Load settings from a configuration file.
        
        Supports JSON format configuration files with validation
        and error handling for malformed or missing files.
        
        Args:
            filepath: Path to the configuration file
            
        Returns:
            True if settings were loaded successfully
            
        Raises:
            ConfigError: If configuration file is invalid or cannot be loaded
        """
        if not os.path.exists(filepath):
            error_msg = f"Configuration file not found: {filepath}"
            logger.error(error_msg)
            raise ConfigError(error_msg, config_key="filepath")
        
        try:
            # Initialize the singleton if not already done
            cls.get_instance()
            
            # Load JSON configuration
            if filepath.endswith('.json'):
                return cls._load_from_json(filepath)
            else:
                error_msg = f"Unsupported config file format: {filepath}"
                logger.error(error_msg)
                raise ConfigError(error_msg, config_key="file_format")
                
        except Exception as e:
            error_msg = f"Error loading settings from {filepath}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def _load_from_json(cls, filepath: str) -> bool:
        """
        Load settings from a JSON file.
        
        Args:
            filepath: Path to the JSON file
            
        Returns:
            True if settings were loaded successfully
        """
        logger.debug(f"Loading configuration from JSON file: {filepath}")
        
        try:
            with open(filepath, 'r') as f:
                settings = json.load(f)
            
            # Update class attributes based on JSON
            for key, value in settings.items():
                if key.startswith('#'):  # Skip comment keys
                    continue
                    
                if hasattr(cls, key):
                    old_value = getattr(cls, key)
                    setattr(cls, key, value)
                    logger.debug(f"Updated config: {key} = {value} (was: {old_value})")
                else:
                    logger.debug(f"Ignoring unknown config key: {key}")
            
            logger.debug(f"Successfully loaded configuration from {filepath}")
            return True
            
        except json.JSONDecodeError as e:
            error_msg = f"Invalid JSON format in {filepath}: {str(e)}"
            logger.error(error_msg)
            raise ConfigError(error_msg, config_key="json_format") from e
        except Exception as e:
            error_msg = f"Error loading JSON settings from {filepath}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def save_to_file(cls, filepath: str) -> bool:
        """
        Save current settings to a JSON file.
        
        Args:
            filepath: Path to save the configuration file
            
        Returns:
            True if settings were saved successfully
            
        Raises:
            ConfigError: If file cannot be written
        """
        logger.debug(f"Saving configuration to file: {filepath}")
        
        try:
            settings = cls.get_all_settings()
            
            with open(filepath, 'w') as f:
                json.dump(settings, f, indent=4)
            
            logger.info(f"Successfully saved configuration to {filepath}")
            return True
            
        except Exception as e:
            error_msg = f"Error saving settings to {filepath}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ConfigError(error_msg) from e
    
    @classmethod
    def get_all_settings(cls) -> Dict[str, Any]:
        """
        Get all settings as a dictionary.
        
        Returns:
            Dictionary of all serializable settings
        """
        settings = {}
        
        # Add all class variables that don't start with underscore
        for key in dir(cls):
            if not key.startswith('_') and not callable(getattr(cls, key)):
                value = getattr(cls, key)
                # Only include serializable types
                if isinstance(value, (str, int, float, bool, list, dict, tuple)) or value is None:
                    settings[key] = value
        
        return settings
    
    @classmethod
    def get_hdbscan_params(cls) -> Dict[str, Any]:
        """
        Get HDBSCAN clustering parameters.
        
        Returns:
            Dictionary of HDBSCAN parameters ready for clustering
        """
        return {
            'min_cluster_size': cls.HDBSCAN_MIN_CLUSTER_SIZE,
            'min_samples': cls.HDBSCAN_MIN_SAMPLES,
            'cluster_selection_epsilon': cls.HDBSCAN_EPSILON,
            'metric': cls.HDBSCAN_METRIC,
            'cluster_selection_method': cls.HDBSCAN_CLUSTER_SELECTION_METHOD,
            'core_dist_n_jobs': 1  # Use single core for reproducibility
        }
    
    @classmethod
    def get_target_tolerance(cls, scale_factor: float = 1.0) -> Dict[str, float]:
        """
        Get target tolerance values with scale factor applied.
        
        Args:
            scale_factor: Scale factor to apply to base tolerance
            
        Returns:
            Dictionary of target names to tolerance values
        """
        # Ensure scale factor is within limits
        scale_factor = max(cls.SCALE_FACTOR_MIN, min(cls.SCALE_FACTOR_MAX, scale_factor))
        
        # Apply scale factor to base tolerance for all targets
        return {target: cls.BASE_TARGET_TOLERANCE * scale_factor 
                for target in cls.EXPECTED_CENTROIDS.keys()}
    
    @classmethod
    def get_plot_dimensions(cls, for_composite: bool = False) -> tuple:
        """
        Get plot dimension settings.
        
        Args:
            for_composite: Whether to get dimensions for composite plot
            
        Returns:
            Figure size as (width, height) tuple
        """
        if for_composite:
            return cls.COMPOSITE_PLOT_SIZE
        else:
            return cls.INDIVIDUAL_FIGURE_SIZE
    
    @classmethod
    def get_axis_limits(cls) -> Dict[str, tuple]:
        """
        Get axis limit settings.
        
        Returns:
            Dictionary of axis limits for x and y axes
        """
        return {
            'x': (cls.X_AXIS_MIN, cls.X_AXIS_MAX),
            'y': (cls.Y_AXIS_MIN, cls.Y_AXIS_MAX)
        }
    
    @classmethod
    def get_grid_intervals(cls) -> Dict[str, int]:
        """
        Get grid interval settings.
        
        Returns:
            Dictionary of grid intervals for x and y axes
        """
        return {
            'x': cls.X_GRID_INTERVAL,
            'y': cls.Y_GRID_INTERVAL
        }

    @classmethod
    def load_user_parameters(cls):
        """
        Load user parameters from the parameter editor if they exist.
        
        This method is called automatically during configuration setup
        to apply any user-customized parameters.
        """
        try:
            from ..utils.parameter_editor import load_parameters_if_exist
            return load_parameters_if_exist(cls)
        except ImportError:
            logger.debug("Parameter editor module not available")
            return False
        except Exception as e:
            logger.debug(f"Could not load user parameters: {e}")
            return False

    # Target colors
    TARGET_COLORS: dict = {}

    @classmethod
    def _get_target_names(cls):
        """
        Return the current list of target names.
        Tries EXPECTED_COPY_NUMBERS first, then EXPECTED_CENTROIDS.
        Always includes special names present in cls.SPECIAL_COLOR_DEFAULTS.
        """
        names = []
        if hasattr(cls, "EXPECTED_COPY_NUMBERS") and isinstance(cls.EXPECTED_COPY_NUMBERS, dict):
            names = list(cls.EXPECTED_COPY_NUMBERS.keys())
        elif hasattr(cls, "EXPECTED_CENTROIDS") and isinstance(cls.EXPECTED_CENTROIDS, dict):
            names = list(cls.EXPECTED_CENTROIDS.keys())

        # ensure specials are present
        for s in getattr(cls, "SPECIAL_COLOR_DEFAULTS", {}).keys():
            if s not in names:
                names.append(s)

        return names

    @classmethod
    def _reconcile_target_colors(cls):
        """
        Make sure TARGET_COLORS matches the current target names.
        Preserves old colors where possible, assigns new ones as needed.
        """
        from itertools import cycle

        old = dict(getattr(cls, "TARGET_COLORS", {}))
        names = cls._get_target_names()
        new = {}

        # keep special defaults stable
        for name, col in getattr(cls, "SPECIAL_COLOR_DEFAULTS", {}).items():
            if name in names:
                new[name] = old.get(name, col)

        # assign colors to all other names
        used = set(new.values())
        base_palette = list(getattr(cls, "DEFAULT_COLOR_PALETTE", []))
        palette = [c for c in base_palette if c not in used]
        color_iter = cycle(base_palette if base_palette else ["#000000"])

        def next_color():
            if palette:
                return palette.pop(0)
            return next(color_iter)

        for name in names:
            if name in new:
                continue
            if name in old and old[name]:
                new[name] = old[name]
            else:
                new[name] = next_color()

        cls.TARGET_COLORS = new

    @classmethod
    def finalize_colors(cls):
        """
        Public method: call after editing targets (e.g., via parameter editor).
        """
        cls._reconcile_target_colors()