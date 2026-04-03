import os
import cv2
import numpy as np
import tifffile
from PIL import Image
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm


def is_meaningful_tile(tile_array, min_tissue_ratio=0.05, color_threshold=30, edge_threshold=0.01):
    """
    Improved pathological tissue filtering function
    :param tile_array: Image data in numpy array format (RGB)
    :param min_tissue_ratio: Minimum tissue area ratio threshold (default 5%)
    :param color_threshold: Color threshold for distinguishing stained regions
    :param edge_threshold: Edge density threshold for detecting tissue structures
    :return: True (contains pathological tissue) or False (does not contain pathological tissue)
    """
    # Convert to grayscale
    gray = cv2.cvtColor(tile_array, cv2.COLOR_RGB2GRAY)

    # Color filtering: detect stained regions
    _, color_mask = cv2.threshold(
        gray, color_threshold, 255, cv2.THRESH_BINARY)

    # Edge detection: detect tissue structures
    edges = cv2.Canny(gray, 50, 150)
    edge_density = np.mean(edges) / 255  # Edge density

    # Calculate stained area ratio
    tissue_ratio = np.sum(color_mask > 0) / color_mask.size

    # Combined criteria
    return (tissue_ratio >= min_tissue_ratio) and (edge_density >= edge_threshold)


def process_tile_wrapper(args):
    """Multi-threaded processing wrapper function"""
    img, i, j, tile_h, tile_w, h, w, output_dir = args
    y_start = i * tile_h
    y_end = (i + 1) * tile_h if i != (h // tile_h - 1) else h
    x_start = j * tile_w
    x_end = (j + 1) * tile_w if j != (w // tile_w - 1) else w

    tile = img[y_start:y_end, x_start:x_end]

    if is_meaningful_tile(tile):
        output_path = os.path.join(output_dir, f"tile_{i:03d}_{j:03d}.png")
        Image.fromarray(tile).save(output_path)
        return True
    return False


def enhanced_split_tif(
    input_path,
    output_dir="enhanced_output",
    grid_size=(10, 10),  # Initial split into 10x10=100 tiles for debugging
    min_tissue_ratio=0.05,  # Minimum tissue area ratio threshold
    color_threshold=30,     # Color threshold
    edge_threshold=0.01,    # Edge density threshold
    max_workers=8           # Number of threads
):
    """Improved pathological slide splitter"""
    os.makedirs(output_dir, exist_ok=True)

    with tifffile.TiffFile(input_path) as tif:
        img = tif.asarray()
        h, w = img.shape[:2]
        rows, cols = grid_size
        tile_h = h // rows
        tile_w = w // cols

        # Generate task list
        tasks = [
            (img, i, j, tile_h, tile_w, h, w, output_dir)
            for i in range(rows)
            for j in range(cols)
        ]

        # Process with progress bar and multi-threading
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(tqdm(
                executor.map(process_tile_wrapper, tasks),
                total=len(tasks),
                desc="Processing Tiles"
            ))

        print(f"\nNumber of valid images saved: {sum(results)}/{len(tasks)}")


# Usage example
enhanced_split_tif(
    input_path="1806935-01.tif",
    output_dir="filtered_tiles",
    grid_size=(10, 10),  # Initial split into 10x10=100 tiles for debugging
    min_tissue_ratio=0.01,  # Minimum tissue area ratio threshold
    color_threshold=10,     # Color threshold
    edge_threshold=0.005,    # Edge density threshold
    max_workers=8           # Number of threads
)
