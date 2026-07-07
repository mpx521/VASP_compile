#!/usr/bin/env python
"""
make_neb_images.py

用ASE的IDPP方法从初态(POSCAR_initial)和末态(POSCAR_final)
插值生成中间images，直接按VASP NEB要求的目录结构输出：

    00/POSCAR   <- 初态（不弛豫，CI-NEB固定端点）
    01/POSCAR   <- 插值中间态1
    02/POSCAR   <- 插值中间态2
    03/POSCAR   <- 插值中间态3
    04/POSCAR   <- 末态（不弛豫，CI-NEB固定端点）

用法：
    # 插值整个体系（基底也会被插值移动）
    python make_neb_images.py POSCAR1 POSCAR2 -n 3

    # 只插值吸附分子，基底原子在所有image里保持初态坐标不变
    # 例如吸附分子是POSCAR里第37-40号原子（按VESTA/POSCAR显示的序号，从1开始数）
    python make_neb_images.py POSCAR1 POSCAR2 -n 3 --moving-atoms 37-40

    # 支持不连续编号，逗号分隔，也可以和区间混用
    python make_neb_images.py POSCAR1 POSCAR2 -n 3 --moving-atoms 37,38,40-42

    # 更省心的方式：按元素类型选择要插值的原子（推荐），比如吸附物种是C/O/H，
    # 基底是Co，直接告诉脚本哪些元素要动，不用去数原子编号
    python make_neb_images.py POSCAR1 POSCAR2 -n 3 --moving-elements C,O,H

    # 如果这次吸附物种只有C和O（没有H），照样写实际有的元素就行
    python make_neb_images.py POSCAR1 POSCAR2 -n 3 --moving-elements C,O

依赖：
    pip install ase --break-system-packages
    (或者你们HPC上通常有现成module，先试 module load ase 或 module load python/ase)

注意：
    - 初态和末态的原子顺序、原子数必须完全一致（同一套原子，只是坐标不同），
      这是NEB插值的硬性要求，如果原子顺序对不上，插值出来的路径没有物理意义
    - 用 --moving-atoms 时，基底坐标全部取自初态(initial)文件，如果你的初态、
      末态文件里基底坐标本身有细微差别（比如两次弛豫收敛程度不同），这里会
      直接忽略末态的基底坐标，统一用初态的。如果这个差别较大，说明你的初末
      态结构本身可能有问题，建议先检查
    - IDPP插值比线性插值更合理，但复杂反应路径（比如涉及键断裂重组、
      或者初末态构型差异很大）时，IDPP插值出来的中间构型仍然可能不合理，
      建议插值完之后用VESTA等软件挨个检查一遍，确认没有原子重叠/键长异常
    - 生成的00和04目录里的POSCAR，在正式提交NEB计算前，建议你分别对初态、
      末态先各自做一次单点或弛豫计算，把对应的OUTCAR拷贝进00/04目录
      （VASP NEB默认要读端点的能量，这一步很多人会漏掉）
"""

import argparse
import os
import sys


def parse_indices(spec):
    """解析形如 '37-40,50,55-56' 的字符串，返回0-indexed的原子序号列表（去重、排序）"""
    indices = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = part.split("-")
            start, end = int(start), int(end)
            indices.update(range(start - 1, end))  # 用户输入按1-indexed
        else:
            indices.add(int(part) - 1)
    return sorted(indices)


def main():
    parser = argparse.ArgumentParser(description="生成CI-NEB中间images")
    parser.add_argument("initial", help="初态结构文件 (如 POSCAR1)")
    parser.add_argument("final", help="末态结构文件 (如 POSCAR2)")
    parser.add_argument("-n", "--nimages", type=int, default=3,
                         help="中间images数量，不含首尾两个端点，默认3")
    parser.add_argument("--linear", action="store_true",
                         help="只用线性插值，不做IDPP优化（默认会做IDPP）")
    parser.add_argument("--moving-atoms", type=str, default=None,
                         help="只插值指定的原子（吸附分子），其余原子（基底）"
                              "在所有image里固定为初态坐标。按POSCAR/VESTA里"
                              "显示的序号(从1开始)填写，支持区间和逗号混合，"
                              "例如 '37-40' 或 '37,38,40-42'。不指定则插值整个体系")
    parser.add_argument("--moving-elements", type=str, default=None,
                         help="按元素符号指定要插值的原子（更省心，推荐），"
                              "逗号分隔，例如 'C,O,H' 或 'C,O'。"
                              "会自动匹配结构里所有该元素的原子，其余元素"
                              "（如基底Co）在所有image里固定为初态坐标。"
                              "和 --moving-atoms 二选一，同时给出时优先用这个")
    args = parser.parse_args()

    try:
        from ase.io import read, write
        from ase.mep import NEB
    except ImportError:
        try:
            from ase.io import read, write
            from ase.neb import NEB
        except ImportError:
            sys.exit(
                "错误：没有找到ASE库。请先安装：\n"
                "  pip install ase --break-system-packages\n"
                "或者在HPC上尝试： module load ase  /  module load python/ase"
            )

    if not os.path.isfile(args.initial):
        sys.exit(f"错误：找不到初态文件 {args.initial}")
    if not os.path.isfile(args.final):
        sys.exit(f"错误：找不到末态文件 {args.final}")

    initial = read(args.initial, format="vasp")
    final = read(args.final, format="vasp")

    if len(initial) != len(final):
        sys.exit(
            f"错误：初态原子数({len(initial)})和末态原子数({len(final)})不一致，"
            "无法插值。请确认两个结构里原子种类、数量、顺序完全对应。"
        )

    n_atoms = len(initial)
    n_total = args.nimages + 2

    if args.moving_elements:
        wanted = set(s.strip() for s in args.moving_elements.split(",") if s.strip())
        symbols = initial.get_chemical_symbols()

        unknown = wanted - set(symbols)
        if unknown:
            sys.exit(f"错误：--moving-elements 里指定的元素 {sorted(unknown)} "
                      f"在结构里找不到。结构里实际有的元素：{sorted(set(symbols))}")

        # 校验初态、末态在这些元素位置上的元素种类是否一致（防止原子顺序错位）
        final_symbols = final.get_chemical_symbols()
        if symbols != final_symbols:
            sys.exit("错误：初态和末态的元素顺序不完全一致，无法按元素筛选插值原子。"
                      "请检查两个POSCAR里原子顺序是否严格对应。")

        moving_idx = [i for i, s in enumerate(symbols) if s in wanted]
        if not moving_idx:
            sys.exit("错误：按指定元素筛选出的原子数为0，请检查元素符号是否写对")

        fixed_elements = sorted(set(symbols) - wanted)
        print(f"体系共 {n_atoms} 个原子，其中元素 {sorted(wanted)} 共 {len(moving_idx)} 个原子"
              f"（吸附分子）将被插值；元素 {fixed_elements} 共 {n_atoms - len(moving_idx)} 个原子"
              f"（基底）在所有image里保持初态坐标不变")

    elif args.moving_atoms:
        moving_idx = parse_indices(args.moving_atoms)
        bad = [i + 1 for i in moving_idx if i < 0 or i >= n_atoms]
        if bad:
            sys.exit(f"错误：--moving-atoms 里有超出范围的序号：{bad}（体系总共{n_atoms}个原子）")

        print(f"体系共 {n_atoms} 个原子，其中 {len(moving_idx)} 个原子（吸附分子）将被插值，"
              f"其余 {n_atoms - len(moving_idx)} 个原子（基底）在所有image里保持初态坐标不变")

    else:
        moving_idx = None

    if moving_idx is not None:
        # 分别对"吸附分子子体系"做插值，基底部分全部锁定为初态坐标
        sub_initial = initial[moving_idx]
        sub_final = final[moving_idx]

        sub_images = [sub_initial]
        sub_images += [sub_initial.copy() for _ in range(args.nimages)]
        sub_images += [sub_final]

        sub_neb = NEB(sub_images)
        sub_neb.interpolate(method="idpp" if not args.linear else "linear")

        # 组装完整体系：基底坐标固定用initial，吸附分子坐标用插值结果覆盖
        images = []
        for i in range(n_total):
            atoms = initial.copy()
            positions = atoms.get_positions()
            positions[moving_idx] = sub_images[i].get_positions()
            atoms.set_positions(positions)
            images.append(atoms)
        # 首尾强制精确等于原始初态/末态（避免子体系插值端点有极小数值误差）
        images[0] = initial.copy()
        images[-1] = final.copy()

    else:
        # 插值整个体系（原有行为）
        images = [initial]
        images += [initial.copy() for _ in range(args.nimages)]
        images += [final]

        neb = NEB(images)
        neb.interpolate(method="idpp" if not args.linear else "linear")

    # 按 00, 01, 02, ..., n_total-1 输出目录
    for i, atoms in enumerate(images):
        dirname = f"{i:02d}"
        os.makedirs(dirname, exist_ok=True)
        outpath = os.path.join(dirname, "POSCAR")
        write(outpath, atoms, format="vasp", direct=True, sort=False)
        tag = "初态(固定)" if i == 0 else ("末态(固定)" if i == n_total - 1 else f"中间image {i}")
        print(f"  已生成 {outpath}   [{tag}]")

    print(f"\n完成。共生成 {n_total} 个目录（含首尾），中间image数：{args.nimages}")
    print("下一步建议：")
    print("  1. 用VESTA等软件挨个打开 01/POSCAR, 02/POSCAR, 03/POSCAR 检查有无原子重叠/异常键长")
    print("  2. 把初态、末态各自单独弛豫/单点计算得到的 OUTCAR 分别拷贝进 00/ 和 04/ 目录")
    print("  3. 在INCAR里设置 IMAGES = {} 并加入NEB相关参数（ICHAIN=0, LCLIMB=.TRUE. 用于CI-NEB等）"
          .format(args.nimages))


if __name__ == "__main__":
    main()
