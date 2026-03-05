# 云端分支与 GitHub `main` 同步：实操清单（中文）

适用场景：
- 你本地已经同步 `main`；
- 云端工作分支（例如 `work`）还有历史提交；
- 合并时提示 notebook 与 `pyotf/microscope.py` 冲突。

---

## 0. 先确认远端与分支状态

```bash
git remote -v
git fetch origin main
git log --oneline --left-right --cherry-pick --no-merges origin/main...HEAD
```

解释：
- 左箭头 `<`：只在 `main` 上；
- 右箭头 `>`：只在当前分支上。

---

## 1. 推荐策略：基于 `main` 进行 rebase

```bash
git rebase origin/main
```

如果出现冲突，Git 会暂停并提示具体文件。

---

## 2. 冲突处理规则（针对你当前项目）

### 2.1 `notebooks/Tutorials/*.ipynb`（强烈建议）

Notebook 是 JSON，容易在 cell id / metadata 上冲突。建议**优先保留 main 版本**：

```bash
git checkout --theirs notebooks/Tutorials/01_项目快速上手与环境检查_zh.ipynb
git checkout --theirs notebooks/Tutorials/02_PSF_OTF_像差与共聚焦针孔测试_zh.ipynb
git checkout --theirs notebooks/Tutorials/03_4Pi与Richards_Wolf参考解工作流_zh.ipynb
```

如需保留当前分支版本，把 `--theirs` 改成 `--ours`。

### 2.2 `pyotf/microscope.py`

代码文件应人工检查冲突块（`<<<<<<<`, `=======`, `>>>>>>>`），按“main 为基线 + 你新增功能”合并。

---

## 3. 标记解决并继续 rebase

```bash
git add notebooks/Tutorials/01_项目快速上手与环境检查_zh.ipynb \
        notebooks/Tutorials/02_PSF_OTF_像差与共聚焦针孔测试_zh.ipynb \
        notebooks/Tutorials/03_4Pi与Richards_Wolf参考解工作流_zh.ipynb \
        pyotf/microscope.py

GIT_EDITOR=true git rebase --continue
```

> 如果你在无交互终端中看到 `EDITOR unset`，使用上面的 `GIT_EDITOR=true`。

---

## 4. 完成后验证

```bash
git status
pytest -q
```

确认：
- 工作区干净；
- 测试通过。

---

## 5. 若你不想 rebase，也可用 merge

```bash
git merge origin/main
```

冲突解决步骤与上面一致。

---

## 6. 何时选择 ours / theirs

- 你要“完全对齐 main”：多数文件用 `--theirs`
- 你要“保留当前分支实验改动”：关键文件用 `--ours`
- notebook 通常优先 `--theirs`，然后在 Jupyter 里手工补充内容

---

## 7. 云端与本地一致性的建议流程

1. 先在本地把改动合入 `main` 并通过测试；
2. 云端分支只做最小增量，频繁 `fetch + rebase`；
3. notebook 变更多时，考虑拆成 `.py` 教程脚本 + notebook 仅展示；
4. 提交前跑 `pytest -q`，避免“能合并但不可运行”。

