
# MatterGen 入门

## 安装与使用

首先克隆 MatterGen 仓库到本地：
```
git clone https://github.com/microsoft/mattergen.git
cd mattergen
```

使用 pip 和 uv 安装项目所需的环境依赖：
```
pip install uv
uv venv .venv --python 3.10
source .venv/bin/activate
uv pip install -e .
```

注意，数据集和模型检查点通过 Git Large File Storage (LFS) 提供。
因此在开始项目前需要确保 LFS 已安装：
```
git lfs --version
```

如果 LFS 未安装，可以使用以下命令安装：
```
sudo apt install git-lfs
git lfs install
```

使用 Git LFS 下载模型检查点：
```
git lfs pull -I checkpoints/
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