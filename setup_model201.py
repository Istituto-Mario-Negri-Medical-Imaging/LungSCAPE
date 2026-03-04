#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import shutil
import subprocess
from pathlib import Path
from urllib.request import urlretrieve


# =========================
# CONFIGURATION (EDIT HERE)

# Set the direct URL to the pretrained nnU-Net v2 model ZIP file.
MODEL_URL = "https://zenodo.org/records/18860421/files/Dataset201_Airways.zip?download=1"

# =========================

MODEL_ZIP_NAME = None

def log(msg: str) -> None:
    print(f"[INFO] {msg}")


def err(msg: str) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)


def which_or_fail(cmd: str) -> str:
    """Ensure a required executable is available in PATH."""
    path = shutil.which(cmd)
    if not path:
        raise RuntimeError(
            f"Command not found: '{cmd}'. "
            "Make sure nnUNetv2 is installed and available in PATH."
        )
    return path


def validate_config() -> None:
    """Validate script configuration before starting."""
    if not MODEL_URL or "YOUR-LINK-HERE" in MODEL_URL:
        raise RuntimeError(
            "MODEL_URL is not configured. Please edit the script and set a valid URL."
        )


def infer_zip_name_from_url(url: str) -> str:
    """Infer the local ZIP filename from the URL."""
    zip_name = Path(url.split("?")[0]).name or "model_pretrained.zip"
    if not zip_name.lower().endswith(".zip"):
        zip_name += ".zip"
    return zip_name


def download_file(url: str, out_path: Path) -> None:
    """Download a file from URL to the target path."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    log(f"Downloading model from: {url}")
    log(f"Saving to: {out_path}")
    urlretrieve(url, str(out_path))
    if not out_path.exists() or out_path.stat().st_size == 0:
        raise RuntimeError("Download failed or file is empty.")
    log(f"Download completed ({out_path.stat().st_size / (1024 * 1024):.1f} MB)")


def create_project_dirs(base_dir: Path) -> dict:
    """Create all project folders, including nnU-Net folders and I/O folders."""
    paths = {
        "downloads": base_dir / "downloads",
        "nnUNet_raw": base_dir / "nnUNet_raw",
        "nnUNet_preprocessed": base_dir / "nnUNet_preprocessed",
        "nnUNet_results": base_dir / "nnUNet_results",
        "INPUT": base_dir / "INPUT",
        "OUTPUT": base_dir / "OUTPUT",
    }
    for name, p in paths.items():
        p.mkdir(parents=True, exist_ok=True)
        log(f"Directory ready: {name} -> {p}")
    return paths


def write_activation_script(paths: dict, out_dir: Path) -> Path:
    """
    Create a shell script (Linux/macOS) or PowerShell script (Windows)
    that sets nnU-Net environment variables for future sessions.
    """
    is_windows = os.name == "nt"
    if is_windows:
        script_path = out_dir / "activate_nnunet_paths.ps1"
        content = f"""$env:nnUNet_raw = "{paths['nnUNet_raw']}"
$env:nnUNet_preprocessed = "{paths['nnUNet_preprocessed']}"
$env:nnUNet_results = "{paths['nnUNet_results']}"

Write-Host "nnU-Net variables set for this PowerShell session:"
Write-Host "  nnUNet_raw=$env:nnUNet_raw"
Write-Host "  nnUNet_preprocessed=$env:nnUNet_preprocessed"
Write-Host "  nnUNet_results=$env:nnUNet_results"
"""
    else:
        script_path = out_dir / "activate_nnunet_paths.sh"
        content = f"""#!/usr/bin/env bash
export nnUNet_raw="{paths['nnUNet_raw']}"
export nnUNet_preprocessed="{paths['nnUNet_preprocessed']}"
export nnUNet_results="{paths['nnUNet_results']}"

echo "nnU-Net variables set for this shell session:"
echo "  nnUNet_raw=$nnUNet_raw"
echo "  nnUNet_preprocessed=$nnUNet_preprocessed"
echo "  nnUNet_results=$nnUNet_results"
"""
    script_path.write_text(content, encoding="utf-8")
    if not is_windows:
        script_path.chmod(0o755)
    log(f"Activation script created: {script_path}")
    return script_path


def run_install_pretrained(zip_path: Path, paths: dict) -> None:
    """Install pretrained model using nnUNetv2_install_pretrained_model_from_zip."""
    nnunet_cmd = which_or_fail("nnUNetv2_install_pretrained_model_from_zip")

    env = os.environ.copy()
    env["nnUNet_raw"] = str(paths["nnUNet_raw"])
    env["nnUNet_preprocessed"] = str(paths["nnUNet_preprocessed"])
    env["nnUNet_results"] = str(paths["nnUNet_results"])

    cmd = [nnunet_cmd, str(zip_path)]
    log("Installing pretrained model...")
    log("Command: " + " ".join(cmd))

    result = subprocess.run(cmd, env=env, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Installation failed (exit code {result.returncode}).")

    log("Installation completed successfully.")


def main():
    parser = argparse.ArgumentParser(
        description="Setup nnU-Net project folder, download pretrained airways model (Dataset201), and install it."
    )
    parser.add_argument(
        "base_dir",
        help="Base directory where ZIP, nnU-Net folders, INPUT and OUTPUT will be created",
    )
    args = parser.parse_args()

    validate_config()

    base_dir = Path(args.base_dir).expanduser().resolve()
    base_dir.mkdir(parents=True, exist_ok=True)

    paths = create_project_dirs(base_dir)

    # Set variables for the current process and child processes
    os.environ["nnUNet_raw"] = str(paths["nnUNet_raw"])
    os.environ["nnUNet_preprocessed"] = str(paths["nnUNet_preprocessed"])
    os.environ["nnUNet_results"] = str(paths["nnUNet_results"])

    zip_name = MODEL_ZIP_NAME if MODEL_ZIP_NAME else infer_zip_name_from_url(MODEL_URL)
    zip_path = paths["downloads"] / zip_name

    try:
        activation_script = write_activation_script(paths, base_dir)

        if not zip_path.exists():
            download_file(MODEL_URL, zip_path)
        else:
            log(f"ZIP already exists, skipping download: {zip_path}")

        run_install_pretrained(zip_path, paths)

        print("\n=== DONE ✅ ===")
        print("Pretrained airways model (Dataset201) installed successfully.")
        print("\nFolders created:")
        for k, v in paths.items():
            print(f"  {k}: {v}")

        print("\nPut your input files into:")
        print(f"  {paths['INPUT']}")
        print("Input filenames must look like: case001_0000.nii.gz, case002_0000.nii.gz, ...")

        print("\nPredictions will be written to:")
        print(f"  {paths['OUTPUT']}")

        print("\nTo use nnU-Net in a new terminal session, load the paths first:")
        if os.name == "nt":
            print(f'  PowerShell:  . "{activation_script}"')
        else:
            print(f'  Bash/Zsh:    source "{activation_script}"')

        print("  Then execute run_model201.py to perform the airways segmentation")

    except Exception as e:
        err(str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
