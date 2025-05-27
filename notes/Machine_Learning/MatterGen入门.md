
# MatterGen 入门

## 安装与使用

首先在 wsl 中启动 nvidia/cuda 的容器：
```
docker run -it --user root --name mattergen-1.0.0 d0117ee15b5f /bin/bash
```

并安装相关依赖：
```
apt update && apt upgrade -y
apt-get update && apt-get upgrade -y
apt-get install -y libcurl4-openssl-dev openssl
apt-get install -y --reinstall git
apt install -y python3.10 python3.10-venv python3.10-dev python3-pip git-lfs
git lfs install
```

然后克隆 MatterGen 仓库到本地：
```
git clone https://github.com/microsoft/mattergen.git
cd mattergen
```

使用 pip 和 uv 安装项目所需的环境依赖：
```
pip install uv
uv venv .venv --python 3.10
source .venv/bin/activate
export UV_HTTP_TIMEOUT=600
uv pip install -e . 
```

使用 Git LFS 下载模型检查点：
```
git lfs pull -I checkpoints/
```

最后将固化全局环境变量：
```
echo 'export PATH="/opt/mattergen/.venv/bin:$PATH"' >> /etc/bash.bashrc
echo 'export VIRTUAL_ENV="/opt/mattergen/.venv"' >> /etc/bash.bashrc
```

然后生成镜像：
```
docker commit mattergen-1.0.0 mattergen:v1.0.0
docker save -o mattergen-1.0.0.tar mattergen:v1.0.0
```

上传 `mattergen-1.0.0.tar` 文件到集群，并修改权限：
```
chmod 777 mattergen-1.0.0.tar 
```

通过 Apptainer 将其转换为 `.sif` 镜像：
```
module load apptainer
apptainer build --sandbox mattergen-sandbox docker-archive:///home/changruiwang-ICME/Software/mattergen/mattergen-1.0.0.tar
apptainer build mattergen.sif mattergen-sandbox
```

## 随机生成

要从预训练的基础模型中采样，可以运行以下命令：
```
export MODEL_NAME=mattergen_base
export RESULTS_PATH=results/
mattergen-generate $RESULTS_PATH --pretrained-name=$MODEL_NAME --batch_size=16 --num_batches=1
```
生成的样本将写入 `$RESULTS_PATH` 目录中。

## 属性生成

要生成具有特定属性的材料，可以使用微调模型。
例如，要从磁密度模型中采样，可以运行以下命令：
```
export MODEL_NAME=dft_mag_density
export RESULTS_PATH="results/$MODEL_NAME/"
mattergen-generate $RESULTS_PATH --pretrained-name=$MODEL_NAME --batch_size=16 --properties_to_condition_on="{'dft_mag_density': 0.15}" --diffusion_guidance_factor=2.0
```
可以根据需要调整 `--diffusion_guidance_factor` 参数，以控制生成样本的多样性和真实性。

## 评估

生成的结构可以使用 MatterSim 机器学习力场进行弛豫，并计算新颖性、唯一性、稳定性等指标：
```
git lfs pull -I data-release/alex-mp/reference_MP2020correction.gz --exclude=""
mattergen-evaluate --structures_path=$RESULTS_PATH --relax=True --structure_matcher='disordered' --save_as="$RESULTS_PATH/metrics.json"
```
评估结果将写入 `$RESULTS_PATH/metrics.json` 文件中，并打印到终端。

## 训练与微调

可以使用以下命令训练 MatterGen 基础模型：
```
mattergen-train data_module=mp_20 ~trainer.logger
```

要微调模型，可以使用以下命令：
```
export PROPERTY=dft_mag_density
mattergen-finetune adapter.pretrained_name=mattergen_base data_module=mp_20 +lightning_module/diffusion_module/model/property_embeddings@adapter.adapter.property_embeddings_adapt.$PROPERTY=$PROPERTY ~trainer.logger data_module.properties=["$PROPERTY"]
```
可以根据需要选择数据集中可用的任何属性进行微调。