# 团队 reference 文件接口

最终实验使用 Testing & Validation 成员确认的数据。当前 `development_reference.npz` 只用于开发验证。

NPZ 使用无 pickle 的数组：

| 键 | 形状/类型 | 要求 |
|---|---|---|
| `t_ref` | `(N,)`, float | 严格递增；起点0、终点40；至少200个正输出时间 |
| `y_ref` | `(N,3)`, float | 每行 `[y1,y2,y3]`；有限；首行为 `[1,0,0]` |
| `metadata` | 可选，标量 JSON 字符串 | 来源、方法、容差、独立核验、有效数字说明 |

推荐时间网格是 `[0]` 加上至少200个从 `1e-8` 到40的对数间隔输出点。求解器内部自适应网格不是这个输出网格。

以下代码只演示将团队**已有**数组导出为约定格式，不计算或替换 reference：

```python
import json
import numpy as np

# t_ref, y_ref 来自团队确认的计算；不要把教材 sanity-check 数字填成轨迹。
metadata = {
    "source": "team_supplied",
    "validation_note": "在这里写实际完成的独立核验、求解器、容差及有效数字依据"
}
np.savez_compressed("data/team_reference.npz",
                    t_ref=t_ref, y_ref=y_ref,
                    metadata=json.dumps(metadata, ensure_ascii=False))
```

`load_reference` 只检查文件、维度、有限性、覆盖区间和初值；不会自动认证科学精度。文件 SHA-256 和 metadata 会进入实验 manifest。

CSV/NPZ 样本输出在 `results/<kind>/runs/`。每份 `_samples.npz` 含共同输出点 `t`、数值解 `Y`、对应 `reference_Y`、是否完成全程、最后一个已接受原生节点 `final_native_t/final_native_Y`。失败 run 只保存可覆盖部分且不外推。

各 run 的 JSON 保存参数与标量统计。较长的 Newton 迭代/残差、自适应接受步/尝试步历史保存在同名 `_history.npz` 中，JSON 的 `history_file` 指向该文件，`history_lengths` 记录数组长度。用 `numpy.load(..., allow_pickle=False)` 读取，避免长数组形成过大的单个 GitHub 文件。

如需全部原生节点：`python code/run_all.py --save-trajectories`。完整轨迹写入 `results/<kind>/trajectories/`，默认不提交 Git。
