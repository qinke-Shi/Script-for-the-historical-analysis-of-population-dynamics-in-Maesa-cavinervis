#!/bin/bash
# One-time environment setup on the server (user-level, no sudo needed).
# 1) switch conda to the Tsinghua mirror (bioconda/bconda-forge are reset by
#    the network otherwise)
# 2) create a python-2.7 env "ys" for the original Yang & Smith scripts
# 3) make sure iqtree2 is available (only if not already installed)

cat > ~/.condarc <<'EOF'
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
custom_channels:
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
envs_dirs:
  - ~/conda_envs
  - ~/.conda/envs
EOF

# python 2.7 env for Yang & Smith scripts (pure standard library, no Biopython needed)
conda create -n ys python=2.7 -y

# iqtree2: check first -- do NOT install if a binary already exists
if ! command -v iqtree2 >/dev/null 2>&1 && ! command -v iqtree >/dev/null 2>&1; then
    echo "iqtree not found; installing via bioconda mirror..."
    conda install -n base -c bioconda iqtree -y
else
    echo "iqtree already installed: $(command -v iqtree2 || command -v iqtree)"
fi

echo "Environment ready. Test:  conda activate ys && python --version"
