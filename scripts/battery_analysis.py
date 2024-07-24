#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用法: python battery_analysis.py

该脚本包含一组用于分析电池材料的方法，通过给定的氧化结构来评估其作为电池的潜力。主要功能包括计算在保持电荷平衡情况下可移除或插入的最大正离子数量，以及计算电池的最大容量。
"""

import math
from collections import defaultdict
from pymatgen.core.periodic_table import Element, Species
from pymatgen.core.structure import Composition
import scipy.constants as const

EV_PER_ATOM_TO_J_PER_MOL = const.e * const.N_A
ELECTRON_TO_AMPERE_HOURS = EV_PER_ATOM_TO_J_PER_MOL / 3600

class BatteryAnalyzer:
    """
    一个用于分析电池材料的类，包含多个方法来评估氧化结构的电池潜力。
    """

    def __init__(self, struc_oxid, cation="Li"):
        """
        初始化方法

        参数:
            struc_oxid: 一个包含氧化态的结构对象
            cation: 正离子的符号或元素对象，必须是带正电荷的，可以是1+/2+/3+等
        """
        for site in struc_oxid:
            if not hasattr(site.specie, "oxi_state"):
                raise ValueError("BatteryAnalyzer 需要结构中指定氧化态！")

        self.struc_oxid = struc_oxid
        self.comp = self.struc_oxid.composition  # 简化后续使用

        if not isinstance(cation, Element):
            self.cation = Element(cation)

        self.cation_charge = self.cation.max_oxidation_state

    @property
    def max_cation_removal(self):
        """
        计算在保持电荷平衡的情况下可以移除的最大正离子数量。

        返回:
            整数，最大可移除的正离子数量，取决于晶胞大小
        """
        oxid_pot = sum(
            [
                (Element(spec.symbol).max_oxidation_state - spec.oxi_state) * self.comp[spec]
                for spec in self.comp
                if is_redox_active_intercalation(Element(spec.symbol))
            ]
        )

        oxid_limit = oxid_pot / self.cation_charge
        num_cation = self.comp[Species(self.cation.symbol, self.cation_charge)]

        return min(oxid_limit, num_cation)

    @property
    def max_cation_insertion(self):
        """
        计算在保持电荷平衡的情况下可以插入的最大正离子数量。

        返回:
            整数，最大可插入的正离子数量，取决于晶胞大小
        """
        lowest_oxid = defaultdict(lambda: 2, {"Cu": 1})  # 仅Cu可以降低到1+
        oxid_list = []

        for spec in self.comp:
            element = Element(spec.symbol)
            if is_redox_active_intercalation(element):
                a = self.comp[spec]
                b = spec.oxi_state
                low = lowest_oxid[spec.symbol]
                c = min(e for e in element.oxidation_states if e >= low)
                t = a * (b - c)
                oxid_list.append(t)

        oxid_pot = sum(oxid_list)
        return oxid_pot / self.cation_charge

    def _get_max_cap_ah(self, num_cations):
        """
        计算插入和移除带电正离子时的最大容量，单位为mAh。

        参数:
            num_cations: 正离子的数量

        返回:
            float，最大容量，单位为mAh
        """
        return num_cations * self.cation_charge * ELECTRON_TO_AMPERE_HOURS

    def get_max_capgrav(self, num_cations):
        """
        计算插入和移除带电正离子时的最大容量，单位为mAh/g。

        参数:
            num_cations: 正离子的数量

        返回:
            float，最大比容量，单位为mAh/g
        """
        weight = self.comp.weight
        return self._get_max_cap_ah(num_cations) / (weight / 1000)

    def get_max_capvol(self, volume=None):
        """
        计算插入和移除带电正离子时的最大容量，单位为mAh/cc。

        参数:
            volume: 用于归一化的体积（默认为初始结构的体积）

        返回:
            float，最大体积容量，单位为mAh/cc
        """
        vol = volume if volume else self.struc_oxid.volume
        return self._get_max_cap_ah(self.max_cation_removal) * 1000 * 1e24 / (vol * const.N_A)

    def get_removals_int_oxid(self):
        """
        返回一组去锂化步骤，使得氧化态为整数。

        返回:
            set，整数正离子移除步长
        """
        oxid_els = [Element(spec.symbol) for spec in self.comp if is_redox_active_intercalation(spec)]
        numa = set()
        for oxid_el in oxid_els:
            numa = numa.union(self._get_int_removals_helper(self.comp.copy(), oxid_el, oxid_els, numa))
        num_cation = self.comp[Species(self.cation.symbol, self.cation_charge)]
        return {num_cation - a for a in numa}

    def _get_int_removals_helper(self, spec_amts_oxi, oxid_el, oxid_els, numa):
        """
        一个辅助方法，用于计算整数氧化态的正离子移除步长。

        参数:
            spec_amts_oxi: 一个包含物种及其数量的字典
            oxid_el: 要氧化的元素
            oxid_els: 可能被氧化的元素列表
            numa: 一个包含正离子数量的集合

        返回:
            set，正离子数量的集合
        """
        oxid_old = min([spec.oxi_state for spec in spec_amts_oxi if spec.symbol == oxid_el.symbol])
        oxid_new = math.floor(oxid_old + 1)
        if oxid_new > oxid_el.max_oxidation_state:
            return numa

        spec_old = Species(oxid_el.symbol, oxid_old)
        spec_new = Species(oxid_el.symbol, oxid_new)
        specamt = spec_amts_oxi[spec_old]
        spec_amts_oxi = {sp: amt for sp, amt in spec_amts_oxi.items() if sp != spec_old}
        spec_amts_oxi[spec_new] = specamt
        spec_amts_oxi = Composition(spec_amts_oxi)

        oxi_noA = sum(
            [spec.oxi_state * spec_amts_oxi[spec] for spec in spec_amts_oxi if spec.symbol not in self.cation.symbol]
        )
        a = max(0, -oxi_noA / self.cation_charge)
        numa = numa.union({a})

        if a == 0:
            return numa
        for ox in oxid_els:
            numa = numa.union(self._get_int_removals_helper(spec_amts_oxi.copy(), ox, oxid_els, numa))
        return numa

def is_redox_active_intercalation(element):
    """
    判断元素是否为氧化还原活性元素，并适合用于嵌入材料。

    参数:
        element: Element对象

    返回:
        bool，元素是否为氧化还原活性元素
    """
    ns = [
        "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Nb", "Mo", "W", 
        "Sb", "Sn", "Bi", "Li", "Na", "K", "Be", "Mg", "Ca", "Al", "Ga", "In", "Tl", "Y"
    ]
    return element.symbol in ns
