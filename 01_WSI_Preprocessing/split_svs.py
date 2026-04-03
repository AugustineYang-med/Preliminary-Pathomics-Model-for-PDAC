import os
import cv2
import numpy as np
from PIL import Image
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import openslide

# =========== Configuration parameters (modify as needed) ============================ #
SVS_PATH = "example.svs"          # .svs file to process
OUTPUT_DIR = "tiles_output"       # Output folder for tile results

GRID_SIZE = (10, 10)              # (rows, cols), i.e. 10x10 = 100 tiles
LEVEL = 0                         # OpenSlide pyramid level: 0=full resolution

MIN_TISSUE_RATIO = 0.01           # Tissue ratio threshold (0~1), lower = more permissive
COLOR_THRESHOLD = 10              # Grayscale threshold (<10 considered stained)
EDGE_THRESHOLD = 0.005           # Edge density threshold (0~1), higher = stricter filtering

MAX_WORKERS = 8                   # Number of concurrent threads: use multiprocessing for CPU-bound tasks
SAVE_FORMAT = "png"               # "png" or "jpeg"
# ================================================ #


# ---------- Check if tile contains valid tissue --------------- #
def is_meaningful_tile(
        tile: np.ndarray,
        min_tissue_ratio: float = 0.05,
        color_threshold: int = 30,
        edge_threshold: float = 0.01) -> bool:

    gray = cv2.cvtColor(tile, cv2.COLOR_RGB2GRAY)

    # Color filtering
    _, color_mask = cv2.threshold(
        gray, color_threshold, 255, cv2.THRESH_BINARY)
    tissue_ratio = np.mean(color_mask > 0)

    # Edge density
    edges = cv2.Canny(gray, 50, 150)
    edge_density = np.mean(edges) / 255.0

    return (tissue_ratio >= min_tissue_ratio) and (edge_density >= edge_threshold)


# ---------- Thread wrapper: open slide within each thread ---------- #
def process_tile(args):
    (svs_path, i, j, tile_h, tile_w, full_h, full_w, level,
     save_dir, thr_params) = args

    slide = openslide.OpenSlide(svs_path)
    y0 = i * tile_h
    x0 = j * tile_w
    h = min(tile_h, full_h - y0)
    w = min(tile_w, full_w - x0)
    region = slide.read_region((x0, y0), level, (w, h))
    tile = np.array(region)[:, :, :3]
    slide.close()

    # Tissue filtering
    if is_meaningful_tile(tile, **thr_params):
        fname = f"tile_{i:03d}_{j:03d}.{SAVE_FORMAT}"
        path = os.path.join(save_dir, fname)
        if SAVE_FORMAT.lower() == "jpeg":
            Image.fromarray(tile).save(path, quality=90, optimize=True)
        else:
            Image.fromarray(tile).save(path)
        return True
    return False


# ---------- Main function ---------- #
def split_svs(
        svs_path: str,
        output_dir: str = "tiles_out",
        grid_size=(10, 10),
        level: int = 0,
        min_tissue_ratio=0.05,
        color_threshold=30,
        edge_threshold=0.01,
        max_workers=8):

    os.makedirs(output_dir, exist_ok=True)
    slide0 = openslide.OpenSlide(svs_path)
    full_w, full_h = slide0.dimensions
    slide0.close()

    rows, cols = grid_size
    tile_h = full_h // rows
    tile_w = full_w // cols

    thr_params = dict(min_tissue_ratio=min_tissue_ratio,
                      color_threshold=color_threshold,
                      edge_threshold=edge_threshold)

    tasks = [(svs_path, i, j, tile_h, tile_w, full_h, full_w,
              level, output_dir, thr_params)
             for i in range(rows) for j in range(cols)]

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(tqdm(executor.map(process_tile, tasks),
                            total=len(tasks), desc="SVS Tiling"))

    kept = sum(results)
    print(f"\n✅ Finished!  Valid tiles retained: {kept}/{len(tasks)}")


# ---------- Entry point ---------- #
if __name__ == "__main__":
    split_svs(
        svs_path=SVS_PATH,
        output_dir=OUTPUT_DIR,
        grid_size=GRID_SIZE,
        level=LEVEL,
        min_tissue_ratio=MIN_TISSUE_RATIO,
        color_threshold=COLOR_THRESHOLD,
        edge_threshold=EDGE_THRESHOLD,
        max_workers=MAX_WORKERS,
    )
