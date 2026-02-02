from exactextract import exact_extract
import exactextract
import zarr
import geopandas as gpd
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple
import logging
import logging.handlers
import sys
from multiprocessing import Queue
import atexit

# Get logger for this module (will be configured by main script)
logger = logging.getLogger(__name__)

# Global queue listener (will be set by setup_logging)
_queue_listener = None


def setup_logging(log_file: str = None, log_level: int = logging.INFO):
    """
    Set up logging to both console and file with multiprocessing support.
    Uses QueueHandler/QueueListener for thread-safe logging from worker processes.
    
    This should be called once at the start of the pipeline.
    
    Parameters:
    -----------
    log_file : str, optional
        Path to log file. If None, only logs to console.
    log_level : int
        Logging level (default: logging.INFO)
        
    Returns:
    --------
    Queue
        Log queue for use in worker processes
    """
    global _queue_listener
    
    # Create root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Remove existing handlers to avoid duplicates
    root_logger.handlers.clear()
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Create handlers that will receive logs from the queue
    handlers = []
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    handlers.append(console_handler)
    
    # File handler (if log_file specified)
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)
    
    # Create queue for multiprocessing-safe logging
    log_queue = Queue(-1)  # No size limit
    
    # Create and start QueueListener
    # This listens to the queue and dispatches to handlers
    _queue_listener = logging.handlers.QueueListener(
        log_queue, *handlers, respect_handler_level=True
    )
    _queue_listener.start()
    
    # Add QueueHandler to root logger
    # Main process logs will go through queue too for consistency
    queue_handler = logging.handlers.QueueHandler(log_queue)
    root_logger.addHandler(queue_handler)
    
    # Register cleanup function
    atexit.register(_cleanup_logging)
    
    # Log initial message
    root_logger.info("Logging initialized with multiprocessing support")
    if log_file:
        root_logger.info(f"Logging to file: {log_file}")
    
    return log_queue


def _cleanup_logging():
    """Clean up logging resources on exit."""
    global _queue_listener
    if _queue_listener is not None:
        _queue_listener.stop()
        _queue_listener = None


def _worker_init(log_queue):
    """
    Initialize worker process. Called once per worker via ProcessPoolExecutor initializer.
    Sets up logging to send all logs to the main process via queue.
    
    Parameters:
    -----------
    log_queue : multiprocessing.Queue
        Queue to send log records to
    """
    # Get root logger for worker process
    root_logger = logging.getLogger()
    root_logger.handlers.clear()
    
    # Add queue handler that sends logs to main process
    queue_handler = logging.handlers.QueueHandler(log_queue)
    root_logger.addHandler(queue_handler)
    root_logger.setLevel(logging.INFO)

def extract_single(
    raster_path: str,
    vectors: gpd.GeoDataFrame,
    operation: list[str] = ['mean'],
    extract_strategy: str = 'feature-sequential',
    include_cols: list[str] = None,
    progress: bool = False
):
    """
    Extract statistics from a single raster for given vector features.
    
    Parameters:
    -----------
    raster_path : str
        Path to the raster file
    vectors : gpd.GeoDataFrame
        GeoDataFrame containing the vector features
    operation : list[str]
        Statistics to extract (e.g., ['mean'], ['mean', 'min', 'max'])
    extract_strategy : str
        Strategy for extraction (default: 'feature-sequential')
    include_cols : list[str]
        Columns from vectors to include in output
    progress : bool
        Show progress bar
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with extracted statistics
    """
    current_extract_df = exact_extract(
        raster_path, 
        vectors, 
        ops=operation, 
        strategy=extract_strategy, 
        include_cols=include_cols,
        progress=progress,
        output='pandas'
    )
    
    return current_extract_df


def _extract_worker(args: Tuple[str, gpd.GeoDataFrame, list, str, list]) -> Tuple[str, pd.DataFrame]:
    """
    Worker function for parallel extraction. Processes a single raster file.
    
    Logging is configured via ProcessPoolExecutor initializer, so it's already set up.
    
    Parameters:
    -----------
    args : tuple
        (raster_path, vectors, operation, extract_strategy, include_cols)
        
    Returns:
    --------
    tuple
        (raster_path, extracted_dataframe)
    """
    raster_path, vectors, operation, extract_strategy, include_cols = args
    
    # Get logger (already configured by worker initializer)
    worker_logger = logging.getLogger(__name__)
    
    try:
        result_df = extract_single(
            raster_path=raster_path,
            vectors=vectors,
            operation=operation,
            extract_strategy=extract_strategy,
            include_cols=include_cols,
            progress=False  # Disable progress for parallel workers
        )
        return (raster_path, result_df)
    except Exception as e:
        worker_logger.error(f"Error processing {raster_path}: {str(e)}")
        return (raster_path, None)


def extract_batch_parallel(
    raster_paths: List[str],
    vectors: gpd.GeoDataFrame,
    operation: list[str] = ['mean'],
    extract_strategy: str = 'feature-sequential',
    include_cols: list[str] = None,
    max_workers: int = None,
    log_queue = None
) -> List[Tuple[str, pd.DataFrame]]:
    """
    Extract statistics from multiple rasters in parallel using concurrent.futures.
    
    Uses ProcessPoolExecutor with initializer to configure logging once per worker process.
    
    Parameters:
    -----------
    raster_paths : List[str]
        List of paths to raster files to process
    vectors : gpd.GeoDataFrame
        GeoDataFrame containing the vector features (shared across all rasters)
    operation : list[str]
        Statistics to extract (e.g., ['mean'])
    extract_strategy : str
        Strategy for extraction (default: 'feature-sequential')
    include_cols : list[str]
        Columns from vectors to include in output
    max_workers : int
        Maximum number of concurrent workers (default: CPU count)
    log_queue : multiprocessing.Queue
        Queue for logging from worker processes
        
    Returns:
    --------
    List[Tuple[str, pd.DataFrame]]
        List of (raster_path, result_dataframe) tuples, in order of completion
    """
    logger.info(f"Starting parallel extraction for {len(raster_paths)} rasters with {max_workers or 'auto'} workers")
    
    # Prepare arguments for each worker
    worker_args = [
        (raster_path, vectors, operation, extract_strategy, include_cols)
        for raster_path in raster_paths
    ]
    
    results = []
    
    # Use initializer to configure logging once per worker process
    with ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=_worker_init,
        initargs=(log_queue,)
    ) as executor:
        # Submit all tasks
        future_to_path = {
            executor.submit(_extract_worker, args): args[0] 
            for args in worker_args
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_path):
            raster_path = future_to_path[future]
            try:
                result = future.result()
                if result[1] is not None:
                    results.append(result)
                    logger.info(f"Completed: {raster_path}")
                else:
                    logger.warning(f"Failed: {raster_path}")
            except Exception as e:
                logger.error(f"Exception for {raster_path}: {str(e)}")
    
    logger.info(f"Batch extraction complete: {len(results)}/{len(raster_paths)} successful")
    
    return results