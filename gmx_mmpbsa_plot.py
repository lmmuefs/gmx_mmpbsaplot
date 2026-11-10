#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gmx_mmpbsa_plot.py — gráficos de MM/PBSA e MM/GBSA a partir de saídas do gmx_MMPBSA (Amber).

O script varre as subpastas, detecta os arquivos de cada sistema e abre um
menu para escolher sistemas, gráficos e opções. Ao final grava um relatório
.txt com tudo o que foi feito e os resultados numéricos.

Arquivos reconhecidos em cada subpasta:
    FINAL_RESULTS_MMPBSA.csv   -> energia vs. tempo e médias ± desvio (gmx_MMPBSA -eo)
    FINAL_DECOMP_MMPBSA.csv    -> decomposição por resíduo (gmx_MMPBSA -deo, &decomp)

Uso:
    python gmx_mmpbsa_plot.py                 # menu interativo, pasta atual
    python gmx_mmpbsa_plot.py /caminho/dados  # outra pasta raiz
    python gmx_mmpbsa_plot.py --help          # opções para uso sem menu

Dependências: numpy, pandas, matplotlib.
"""
from __future__ import annotations

import argparse
import configparser
import csv
import datetime as dt
import logging
import re
import sys
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# ============================================================
# CONSTANTES
# ============================================================

KJ_TO_KCAL = 4.184
OUT_DIR_NAME = "MM_PBSA_plots"
CONFIG_NAME = "plot_style.ini"

XVG_NAME = "FINAL_RESULTS_MMPBSA.csv"
RES_NAME = "FINAL_DECOMP_MMPBSA.csv"
RES_DAT_NAME = "FINAL_DECOMP_MMPBSA.dat"   # alternativa: média ± DP já agregados
SUM_NAME = "FINAL_RESULTS_MMPBSA.csv"

# Paleta padrão dos sistemas (distinguível, amigável a daltonismo e à impressão)
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#7F7F7F", "#000000"]

# Cores padrão dos componentes no gráfico por resíduo
COMPONENT_COLORS = ["#18c5d9", "#e600a9", "#9b00d4", "#38e600", "#000000"]

# Componentes: chave interna -> (coluna média, coluna desvio) no residues CSV
COMPONENTS = ["vdw", "electrostatic", "polar", "nonpolar", "total"]
RES_COLS = {
    "vdw": ("vDW", "vdW-stddev"),
    "electrostatic": ("Electrostatic", "Electrostatic-stddev"),
    "polar": ("polar", "polar-stddev"),
    "nonpolar": ("apolar", "apolar-stddev"),
    "total": ("total", "total-stddev"),
}
# Resíduos que NÃO são considerados ligante na detecção automática
STANDARD_RESIDUES = set("""
ALA ARG ASN ASP CYS GLN GLU GLY HIS HID HIE HIP HSD HSE HSP ILE LEU LYS MET PHE
PRO SER THR TRP TYR VAL CYX CYM ASH GLH LYN ACE NME NMA NH2 NA CL K MG ZN CA
SOL WAT HOH
""".split())

# Nomes das pastas de saída, por idioma
GROUP_DIRS = {
    "en": {"time": "energy_vs_time", "summary": "summary",
           "residues": "per_residue", "tables": "tables"},
    "pt": {"time": "energia_vs_tempo", "summary": "medias",
           "residues": "por_residuo", "tables": "tabelas"},
}

# Textos do relatório .txt e das tabelas, por idioma
REPORT_TEXT = {
    "en": dict(
        title="MM/PBSA REPORT", date="Date", root="Root folder",
        out="Output folder", cfg="Configuration", settings="SETTINGS",
        systems="Systems and labels", color="color", plots="Plots",
        window="Time range", residues="Residues", all="All",
        ligand_kept="Ligand kept in ranking", yes="yes", no="no",
        ligand_det="ligand detected in", none="none", formats="Formats",
        layout="Folder layout", lang="Language", font="Font",
        font_sub="substituting",
        layout_names={"type": "by plot type", "format": "by file format",
                      "flat": "single folder"},
        units="Units: energies in kcal/mol (as reported by gmx_MMPBSA); time in ns.",
        time_sec="ENERGY VS. TIME — mean ± SD between {a:g} and {b:g} ns (kcal/mol)",
        system="System", frames="frames",
        avg_sec="COMPONENT AVERAGES — full trajectory, from FINAL_RESULTS_MMPBSA.csv "
                "(kcal/mol; ± = SD over frames)",
        res_sec="PER-RESIDUE DECOMPOSITION (kcal/mol; more negative = more "
                "favorable contribution)",
        full_saved="complete ranking saved as CSV; first 10 below",
        res_window="residues computed only from the frames in {a:g}–{b:g} ns",
        ligand_excl="ligand {l} excluded", ligand_keep="ligand kept",
        files="GENERATED FILES", warnings="WARNINGS",
        sd="± SD", rank="Rank", residue="Residue", complex="Complex",
        pm_cols={"vdw": "vDW ± std", "electrostatic": "Elec. ± std",
                 "polar": "Polar ± std", "nonpolar": "Apolar ± std",
                 "total": "Total ± std"},
    ),
    "pt": dict(
        title="RELATÓRIO MM/PBSA", date="Data", root="Pasta raiz",
        out="Pasta de saída", cfg="Configuração", settings="ESCOLHAS",
        systems="Sistemas e legendas", color="cor", plots="Gráficos",
        window="Faixa de tempo", residues="Resíduos", all="Todos",
        ligand_kept="Ligante mantido no ranking", yes="sim", no="não",
        ligand_det="ligante detectado em", none="nenhum", formats="Formatos",
        layout="Organização das pastas", lang="Idioma", font="Fonte",
        font_sub="substituindo",
        layout_names={"type": "por tipo de gráfico", "format": "por formato",
                      "flat": "tudo na mesma pasta"},
        units="Unidades: energias em kcal/mol (como o gmx_MMPBSA reporta); tempo em ns.",
        time_sec="ENERGIA NO TEMPO — média ± desvio entre {a:g} e {b:g} ns (kcal/mol)",
        system="Sistema", frames="frames",
        avg_sec="MÉDIAS DOS COMPONENTES — trajetória completa, de "
                "FINAL_RESULTS_MMPBSA.csv (kcal/mol; ± = desvio padrão entre frames)",
        res_sec="DECOMPOSIÇÃO POR RESÍDUO (kcal/mol; mais negativo = "
                "contribuição mais favorável)",
        full_saved="ranking completo salvo em CSV; abaixo, os 10 primeiros",
        res_window="resíduos calculados só com os frames em {a:g}–{b:g} ns",
        ligand_excl="ligante {l} excluído", ligand_keep="ligante mantido",
        files="ARQUIVOS GERADOS", warnings="AVISOS",
        sd="± DP", rank="Pos.", residue="Resíduo", complex="Complexo",
        pm_cols={"vdw": "vDW ± dp", "electrostatic": "Elet. ± dp",
                 "polar": "Polar ± dp", "nonpolar": "Apolar ± dp",
                 "total": "Total ± dp"},
    ),
}

TEXT_EN = {
    "time_xlabel": "TIME (ns)",
    "energy_ylabel": "ENERGY (kcal/mol)",
    "residue_xlabel": "RESIDUE",
    "residue_number_xlabel": "RESIDUE NUMBER",
    "residue_title": "Per-residue energy decomposition",
    "profile_title": "Per-residue total energy profile",
    "summary_title": "Average energy components (full trajectory)",
    "title_total": "Total Binding Energy",
    "title_vdw": "Van der Waals",
    "title_electrostatic": "Electrostatics",
    "title_polar": "Polar Solvation",
    "title_nonpolar": "Nonpolar Solvation",
    "comp_vdw": "vDW",
    "comp_electrostatic": "Electrostatic",
    "comp_polar": "Polar",
    "comp_nonpolar": "Nonpolar",
    "comp_total": "Total",
    "title_gas": "Gas-Phase Energy (ΔG gas)",
    "title_solv": "Solvation Energy (ΔG solv)",
    "comp_gas": "Gas (GGAS)",
    "comp_solv": "Solvation (GSOLV)",
}
TEXT_PT = {
    "time_xlabel": "TEMPO (ns)",
    "energy_ylabel": "ENERGIA (kcal/mol)",
    "residue_xlabel": "RESÍDUO",
    "residue_number_xlabel": "NÚMERO DO RESÍDUO",
    "residue_title": "Decomposição de energia por resíduo",
    "profile_title": "Perfil de energia total por resíduo",
    "summary_title": "Média dos componentes de energia (trajetória completa)",
    "title_total": "Energia de Ligação Total",
    "title_vdw": "van der Waals",
    "title_electrostatic": "Eletrostática",
    "title_polar": "Solvatação Polar",
    "title_nonpolar": "Solvatação Apolar",
    "comp_vdw": "vDW",
    "comp_electrostatic": "Eletrostática",
    "comp_polar": "Polar",
    "comp_nonpolar": "Apolar",
    "comp_total": "Total",
    "title_gas": "Energia em Fase Gasosa (ΔG gas)",
    "title_solv": "Energia de Solvatação (ΔG solv)",
    "comp_gas": "Gás (GGAS)",
    "comp_solv": "Solvatação (GSOLV)",
}

# Configuração editável (plot_style.ini)
DEFAULT_CONFIG: Dict[str, Dict[str, str]] = {
    "general": {"language": "en", "version": "2"},
    "names": {},
    "text_en": TEXT_EN,
    "text_pt": TEXT_PT,
    "font": {
        "family": "Arial",
        "size": "12",
        "title_size": "14",
        "label_size": "12",
        "tick_size": "10",
        "legend_size": "10",
    },
    "style": {
        "dpi": "600",
        "line_width": "1.3",
    },
}

CONFIG_HEADER = """\
# Configuração do gmx_mmpbsa_plot.py — edite à vontade e rode o script de novo.
#
# [general]  : language = en ou pt (último idioma usado; vale para o menu,
#              gráficos, tabelas e relatório).
# [names]    : nome de legenda de cada sistema (chave = pasta). Pode usar
#              mathtext do matplotlib, ex.:  holo1 = $\\it{Mtb}$PS - KS1
# [text_en] / [text_pt] : títulos, rótulos de eixo e nomes dos componentes,
#              um bloco por idioma.
# [font]     : fonte e tamanhos.   [style] : dpi e espessura de linha.
#
# Obs.: ao regravar, o script reescreve este arquivo e apaga comentários
# adicionais que você inserir; este cabeçalho é sempre recolocado.

"""


UI_LANG = ["pt"]      # idioma da interface (menus, avisos e mensagens)


def _(msg: str, /, **kw) -> str:
    """Traduz uma mensagem da interface (as chaves são os textos em português)."""
    out = EN_UI.get(msg, msg) if UI_LANG[0] == "en" else msg
    return out.format(**kw) if kw else out


def default_ui_lang() -> str:
    """Idioma sugerido na primeira vez: o do sistema (pt* -> pt, resto -> en)."""
    import os
    loc = (os.environ.get("LC_ALL") or os.environ.get("LC_MESSAGES")
           or os.environ.get("LANG") or "")
    return "pt" if loc.lower().startswith("pt") else "en"


class DataError(Exception):
    """Arquivo de dados em formato inesperado."""


class UsageError(Exception):
    """Pedido inconsistente do usuário."""


# ============================================================
# ESTRUTURAS
# ============================================================

@dataclass
class System:
    key: str                     # caminho relativo da pasta (identificador)
    path: Path
    xvg: Optional[Path] = None
    res: Optional[Path] = None
    summary: Optional[Path] = None
    label: str = ""
    color: str = ""
    ligand: Optional[str] = None
    drop_ligand: Optional[str] = None

    @property
    def slug(self) -> str:
        return re.sub(r"[^\w.-]+", "_", self.key).strip("_") or "root"


class Report:
    """Acumula o texto do relatório final."""

    def __init__(self) -> None:
        self.lines: List[str] = []
        self.warnings: List[str] = []
        self.files: List[Path] = []
        self.base: Optional[Path] = None      # para exibir caminhos curtos

    def add(self, text: str = "") -> None:
        self.lines.append(text)

    def section(self, title: str) -> None:
        self.add()
        self.add("=" * 78)
        self.add(title)
        self.add("=" * 78)

    def warn(self, text: str) -> None:
        print(f"  ⚠ {text}")
        self.warnings.append(text)

    def file(self, path: Path) -> None:
        self.files.append(path)
        try:
            shown = path.relative_to(self.base) if self.base else path
        except ValueError:
            shown = path
        print(f"    ✓ {shown}")


@dataclass
class Ctx:
    root: Path
    outdir: Path
    cfg: configparser.ConfigParser
    report: Report
    formats: List[str] = field(default_factory=lambda: ["png", "pdf", "svg"])
    interactive: bool = False
    export_csv: bool = False
    lang: str = "en"
    layout: str = "type"          # type | format | flat

    def t(self, key: str) -> str:
        try:
            return self.cfg[f"text_{self.lang}"][key]
        except KeyError:
            return (TEXT_PT if self.lang == "pt" else TEXT_EN)[key]

    def target_dir(self, group: str, fmt: Optional[str] = None) -> Path:
        """Pasta de destino conforme a organização escolhida."""
        names = GROUP_DIRS[self.lang]
        if group == "tables":
            d = self.outdir / names["tables"]
        elif self.layout == "type":
            d = self.outdir / names[group]
        elif self.layout == "format" and fmt:
            d = self.outdir / fmt
        else:
            d = self.outdir
        d.mkdir(parents=True, exist_ok=True)
        return d


# ============================================================
# ENTRADA DO USUÁRIO
# ============================================================

AUTO = False   # --yes: usa padrões e nunca pergunta
_STEP = [0]    # contador das etapas exibidas no menu

USE_COLOR = sys.stdout.isatty()


def _style(code: str, text: str) -> str:
    return f"\033[{code}m{text}\033[0m" if USE_COLOR else text


def bold(text: str) -> str:
    return _style("1", text)


def dim(text: str) -> str:
    return _style("2", text)


def step(title: str) -> None:
    """Cabeçalho numerado de cada etapa do menu (omitido no modo --yes)."""
    if AUTO:
        return
    _STEP[0] += 1
    print()
    print(bold(f"{_STEP[0]} · {_(title)}"))
    print(dim("─" * 64))


def ask(text: str, parse, preset: Optional[str] = None,
        default: Optional[str] = None):
    """Pergunta até receber algo que `parse` aceite (ValueError = inválido)."""
    text = _(text)
    if preset is not None:
        try:
            return parse(preset)
        except ValueError as e:
            raise UsageError(_("Argumento inválido para '{text}': {err}",
                               text=text, err=e))
    if AUTO:
        if default is None:
            raise UsageError(_("Modo --yes: '{text}' não tem padrão; "
                               "informe por argumento.", text=text))
        return parse(default)
    shown = f" [{default}]" if default not in (None, "") else ""
    while True:
        raw = input(f"  {text}{shown}: ").strip()
        if not raw and default is not None:
            raw = default
        try:
            return parse(raw)
        except ValueError as e:
            print(f"  ! {e}")


def parse_yn(s: str) -> bool:
    s = s.strip().lower()
    if s in ("s", "sim", "y", "yes"):
        return True
    if s in ("n", "nao", "não", "no", ""):
        return False
    raise ValueError(_("responda s ou n"))


def parse_indices(text: str, n: int) -> List[int]:
    text = text.strip().lower()
    if text in ("all", "todos", "a", "*"):
        return list(range(n))
    out: List[int] = []
    for part in re.split(r"[\s,;]+", text):
        if not part:
            continue
        if re.fullmatch(r"\d+-\d+", part):
            a, b = (int(x) for x in part.split("-"))
            rng = range(a, b + 1)
        elif part.isdigit():
            rng = [int(part)]
        else:
            raise ValueError(_("'{part}' não é um número válido", part=part))
        for k in rng:
            if not 1 <= k <= n:
                raise ValueError(_("{k} fora da faixa 1-{n}", k=k, n=n))
            if k - 1 not in out:
                out.append(k - 1)
    if not out:
        raise ValueError(_("nada selecionado"))
    return out


def parse_hex_list(text: str) -> List[str]:
    out = []
    for part in re.split(r"[\s,;]+", text.strip()):
        if not part:
            continue
        if not re.fullmatch(r"#?[0-9a-fA-F]{6}", part):
            raise ValueError(_("'{part}' não é um hex válido (ex.: #1f77b4)", part=part))
        out.append("#" + part.lstrip("#").lower())
    return out


def parse_formats(text: str) -> List[str]:
    text = text.strip().lower()
    if not text or text in ("all", "todos"):
        return ["png", "pdf", "svg"]
    out = []
    for part in re.split(r"[\s,;]+", text):
        if part not in ("png", "pdf", "svg"):
            raise ValueError(_("formato '{part}' não suportado (png, pdf, svg)", part=part))
        if part not in out:
            out.append(part)
    return out


def parse_pos_float(text: str) -> float:
    try:
        v = float(text.replace(",", "."))
    except ValueError:
        raise ValueError(_("informe um número, ex.: 0.01"))
    if v <= 0:
        raise ValueError(_("deve ser maior que zero"))
    return v


def parse_window(text: str) -> Tuple[float, float]:
    nums = re.findall(r"\d+(?:\.\d+)?", text)
    if len(nums) != 2:
        raise ValueError(_("informe dois números: início e fim (ex.: 20 100)"))
    a, b = float(nums[0]), float(nums[1])
    if a >= b:
        raise ValueError(_("o início precisa ser menor que o fim"))
    return a, b


def parse_top(text: str) -> Optional[int]:
    """Retorna 5, 10 ou None (= todos)."""
    s = text.strip().lower()
    if s in ("1", "5", "top5"):
        return 5
    if s in ("2", "10", "top10"):
        return 10
    if s in ("3", "all", "todos"):
        return None
    raise ValueError(_("escolha 1 (Top 5), 2 (Top 10) ou 3 (Todos)"))


def parse_top_limited(text: str) -> int:
    n = parse_top(text)
    if n is None:
        raise ValueError(_("escolha 1 (Top 5) ou 2 (Top 10)"))
    return n


def parse_lang(text: str) -> str:
    t = text.strip().lower()
    if t in ("1", "pt", "pt-br", "ptbr", "português", "portugues"):
        return "pt"
    if t in ("2", "en", "en-us", "enus", "eng", "english", "inglês", "ingles"):
        return "en"
    raise ValueError("escolha 1 (Português) ou 2 (English) / choose 1 or 2")


def parse_layout(text: str) -> str:
    t = text.strip().lower()
    if t in ("1", "type", "tipo"):
        return "type"
    if t in ("2", "format", "formato"):
        return "format"
    if t in ("3", "flat", "junto"):
        return "flat"
    raise ValueError(_("escolha 1 (por tipo), 2 (por formato) ou 3 (tudo junto)"))


def parse_name(text: str) -> str:
    name = text.strip()
    if not name:
        raise ValueError(_("nome vazio"))
    return name


# ============================================================
# CONFIGURAÇÃO (plot_style.ini)
# ============================================================

def peek_saved_lang(path: Path) -> str:
    """Último idioma usado (lido do .ini) ou o do sistema."""
    try:
        cp = configparser.ConfigParser(interpolation=None)
        cp.read(path, encoding="utf-8")
        v = cp.get("general", "language", fallback="").strip().lower()
        if v in ("pt", "en"):
            return v
    except (configparser.Error, OSError):
        pass
    return default_ui_lang()


def load_config(path: Path) -> Tuple[configparser.ConfigParser, bool]:
    cfg = configparser.ConfigParser(interpolation=None)
    cfg.optionxform = str          # preserva maiúsculas/minúsculas das chaves
    changed = False
    if path.exists():
        try:
            cfg.read(path, encoding="utf-8")
        except configparser.Error as e:
            raise DataError(_("{path} está malformado: {err}", path=path, err=e))
    else:
        changed = True
    legacy = path.exists() and not (cfg.has_section("general")
                                    and cfg.has_option("general", "version"))
    if legacy:                       # formato antigo: bloco único [text], fonte serif
        if cfg.has_section("text"):
            if not cfg.has_section("text_en"):
                cfg.add_section("text_en")
            for k, v in cfg["text"].items():
                if k in TEXT_EN and not cfg.has_option("text_en", k):
                    cfg.set("text_en", k, v)
            cfg.remove_section("text")
        if cfg.has_section("font") and cfg["font"].get("family") == "Times New Roman":
            cfg["font"]["family"] = "Arial"
        if not cfg.has_section("general"):
            cfg.add_section("general")
        cfg.set("general", "version", "2")
        changed = True
        print(dim(_("  Configuração antiga atualizada (fonte Arial, textos por idioma).")))
    for sec, items in DEFAULT_CONFIG.items():
        if not cfg.has_section(sec):
            cfg.add_section(sec)
            changed = True
        for k, v in items.items():
            if not cfg.has_option(sec, k):
                cfg.set(sec, k, v)
                changed = True
    return cfg, changed


def save_config(cfg: configparser.ConfigParser, path: Path) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(CONFIG_HEADER)
        cfg.write(fh)


def cfg_float(cfg: configparser.ConfigParser, sec: str, key: str) -> float:
    try:
        return float(cfg[sec][key])
    except (ValueError, KeyError):
        return float(DEFAULT_CONFIG[sec][key])


def resolve_font(family: str) -> str:
    """Nome da fonte realmente usada (família pedida ou a primeira reserva)."""
    for name in (family, "Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"):
        try:
            path = font_manager.findfont(name, fallback_to_default=False)
            return font_manager.FontProperties(fname=path).get_name()
        except ValueError:
            continue
    return "DejaVu Sans"


def apply_style(cfg: configparser.ConfigParser) -> None:
    fam = cfg["font"].get("family", "Arial")
    plt.rcParams.update({
        "font.family": "sans-serif",
        # Liberation Sans tem as mesmas métricas do Arial (reserva em Linux)
        "font.sans-serif": [fam, "Arial", "Liberation Sans", "Helvetica",
                            "DejaVu Sans"],
        # mathtext (ex.: $\it{Mtb}$) na mesma família da fonte do gráfico
        "mathtext.fontset": "custom",
        "mathtext.rm": "sans",
        "mathtext.it": "sans:italic",
        "mathtext.bf": "sans:bold",
        "mathtext.default": "regular",
        "font.size": cfg_float(cfg, "font", "size"),
        "axes.titlesize": cfg_float(cfg, "font", "title_size"),
        "axes.labelsize": cfg_float(cfg, "font", "label_size"),
        "xtick.labelsize": cfg_float(cfg, "font", "tick_size"),
        "ytick.labelsize": cfg_float(cfg, "font", "tick_size"),
        "legend.fontsize": cfg_float(cfg, "font", "legend_size"),
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        # texto continua texto nos arquivos vetoriais (editável no Inkscape etc.)
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


# ============================================================
# DESCOBERTA E LEITURA DE ARQUIVOS
# ============================================================

def discover(root: Path, outdir: Path) -> List[System]:
    found: Dict[Path, System] = {}
    wanted = {XVG_NAME: "xvg", RES_NAME: "res", RES_DAT_NAME: "res"}
    for f in sorted(root.rglob("*")):
        if f.name not in wanted or not f.is_file():
            continue
        rel_parts = f.relative_to(root).parts
        if any(p.startswith(".") for p in rel_parts):
            continue
        if outdir in f.parents or any(p.startswith(OUT_DIR_NAME) for p in rel_parts):
            continue
        d = f.parent
        if d not in found:
            rel = d.relative_to(root).as_posix()
            found[d] = System(key=rel if rel != "." else d.name, path=d)
        if f.name == RES_DAT_NAME and found[d].res is not None:
            continue                      # o CSV por frame tem prioridade
        setattr(found[d], wanted[f.name], f)
        if f.name == XVG_NAME:
            found[d].summary = f
    return sorted(found.values(), key=lambda s: _natural_key(s.key))


def _natural_key(s: str):
    return [int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", s)]


def _num(x: str) -> Optional[float]:
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def parse_energy_csv(path: Path, model: Optional[str]) -> Optional[pd.DataFrame]:
    """Lê a tabela 'Delta Energy Terms' (por frame) do CSV do gmx_MMPBSA."""
    rows = list(csv.reader(path.read_text(errors="replace").splitlines()))
    tables: Dict[str, List[List[str]]] = {}
    headers: Dict[str, List[str]] = {}
    cur_model, in_delta, header = "", False, None
    for r in rows:
        first = r[0].strip() if r else ""
        low = first.lower()
        if not first:
            if in_delta and header is not None and tables.get(cur_model):
                in_delta = False
            continue
        if "generalized born" in low or low in ("gb", "|gb"):
            cur_model = "GB"
        elif "poisson" in low or low in ("pb", "|pb"):
            cur_model = "PB"
        if re.search(r"delta.*energy terms|^deltas?\b", low) and not _num(first):
            in_delta, header = True, None
            continue
        if re.search(r"(complex|receptor|ligand).*energy terms", low):
            in_delta = False
            continue
        if not in_delta:
            continue
        if header is None:
            if low.startswith("frame"):
                header = [c.strip().upper() for c in r]
                headers[cur_model] = header
                tables.setdefault(cur_model, [])
            continue
        if _num(first) is not None:
            tables[cur_model].append(r)
    if not tables:
        return None
    if model and model in tables:
        key = model
    else:
        key = next(iter(tables))
    if key:
        USED_MODELS.add(key)
    df = pd.DataFrame(tables[key], columns=headers[key][: len(tables[key][0])]).apply(pd.to_numeric, errors="coerce")
    col = {c.upper(): c for c in df.columns}

    def pick(*names):
        for n in names:
            if n in col:
                return df[col[n]]
        return None

    out = pd.DataFrame({"frame": df[col[next(iter(col))]].astype(int)})
    vdw, eel = pick("VDWAALS"), pick("EEL")
    pol = pick("EPB", "EGB", "EPB(NPT)")
    gsolv, total, gas = pick("GSOLV"), pick("TOTAL"), pick("GGAS")
    if vdw is None or eel is None or pol is None or total is None:
        return None
    out["VDW"], out["EEL"], out["POLAR"] = vdw, eel, pol
    # apolar: GB -> ESURF; PB -> ENPOLAR (+ EDISPER) ou ECAVITY + EDISPER; reserva: GSOLV - polar
    if pick("ESURF") is not None:
        nonpolar = pick("ESURF")
    elif pick("ENPOLAR") is not None:
        nonpolar = pick("ENPOLAR") + (pick("EDISPER") if pick("EDISPER") is not None else 0.0)
    elif pick("ECAVITY") is not None:
        nonpolar = pick("ECAVITY") + (pick("EDISPER") if pick("EDISPER") is not None else 0.0)
    elif gsolv is not None:
        nonpolar = gsolv - pol
    else:
        return None
    out["NONPOLAR"] = nonpolar.round(4)
    out["GAS"] = gas if gas is not None else (vdw + eel)
    out["SOLV"] = gsolv if gsolv is not None else (pol + nonpolar)
    out["TOTAL"] = total
    return out.reset_index(drop=True)


def parse_decomp_csv(path: Path) -> Optional[pd.DataFrame]:
    """Lê 'Total Energy Decomposition' da seção DELTAS do FINAL_DECOMP_MMPBSA.csv."""
    rows = list(csv.reader(path.read_text(errors="replace").splitlines()))
    in_delta, in_total, header = False, False, None
    data: List[List[str]] = []
    for r in rows:
        first = r[0].strip() if r else ""
        low = first.lower()
        if not first and not (len(r) > 1 and any(c.strip() for c in r)):
            if in_total and data:
                break
            continue
        if low.startswith("delta"):
            in_delta, in_total, header = True, False, None
            continue
        if low.startswith(("complex", "receptor", "ligand")) and not _num(first):
            in_delta = False
            continue
        if not in_delta:
            continue
        if "total energy decomposition" in low:
            in_total, header = True, None
            continue
        if re.search(r"(sidechain|backbone).*decomposition", low):
            in_total = False
            continue
        if not in_total:
            continue
        if header is None:
            if low.startswith("frame") or "residue" in [c.strip().lower() for c in r]:
                header = [c.strip().lower() for c in r]
            continue
        if len(r) >= len(header) and _num(first) is not None:
            data.append(r[: len(header)])
    if not data or header is None:
        return None
    df = pd.DataFrame(data, columns=header)

    def find(*keys):
        for c in df.columns:
            if any(k in c for k in keys):
                return c
        return None

    res_c = find("residue")
    cols = {
        "INTERNAL": find("internal"), "VDW": find("van der waals"), "EEL": find("electrostatic"),
        "POLAR": find("polar solv"), "NONPOLAR": find("non-polar", "nonpolar"), "TOTAL": find("total"),
    }
    if res_c is None or any(v is None for v in cols.values()):
        return None
    out = pd.DataFrame({"frame": pd.to_numeric(df[df.columns[0]], errors="coerce"), "residue": df[res_c].str.strip()})
    for k, c in cols.items():
        out[k] = pd.to_numeric(df[c], errors="coerce")
    parts = out["residue"].str.split(":", expand=True)
    # formato gmx_MMPBSA: R:A:LEU:25  (R/L = receptor/ligante)
    if parts.shape[1] >= 4:
        out["mol"], out["chain"], out["resname"], out["resnum"] = parts[0], parts[1], parts[2], pd.to_numeric(parts[3], errors="coerce")
    else:
        out["mol"] = ""
        out["chain"] = ""
        out["resname"] = out["residue"].str.extract(r"([A-Za-z]+)")[0]
        out["resnum"] = pd.to_numeric(out["residue"].str.extract(r"(\d+)")[0], errors="coerce")
    out["label"] = out["resname"].str.upper() + out["resnum"].astype("Int64").astype(str)
    return out.dropna(subset=["resnum"]).reset_index(drop=True)



def parse_decomp_dat(path: Path) -> Optional[pd.DataFrame]:
    """FINAL_DECOMP_MMPBSA.dat: 'Total Energy Decomposition' da seção DELTAS (média, DP, EPM)."""
    in_delta = in_total = False
    rows = []
    for ln in path.read_text(errors="replace").splitlines():
        t = ln.strip()
        low = t.lower()
        if low.startswith("deltas"):
            in_delta, in_total = True, False
        elif low.startswith(("complex", "receptor", "ligand")) and low.endswith(":"):
            in_delta = in_total = False
        elif in_delta and "total energy decomposition" in low:
            in_total = True
        elif in_delta and re.search(r"(sidechain|backbone) energy decomposition", low):
            in_total = False
        elif in_total and not t:
            if rows:
                break
        elif in_total:
            f = t.split(",")
            if len(f) >= 19 and f[0].count(":") >= 3 and _num(f[1]) is not None:
                rows.append(f)
    if not rows:
        return None
    names = ["INTERNAL", "VDW", "EEL", "POLAR", "NONPOLAR", "TOTAL"]
    out = []
    for f in rows:
        mol, chain, resname, resnum = f[0].split(":")[:4]
        d = {"mol": mol, "chain": chain, "resname": resname, "resnum": _num(resnum)}
        for i, n in enumerate(names):
            d[f"m_{n}"], d[f"sd_{n}"] = float(f[1 + 3 * i]), float(f[2 + 3 * i])
        out.append(d)
    return pd.DataFrame(out)


NS_PER_FRAME = [None]   # ns por frame da trajetória (coluna 'Frame #' do CSV)
MODEL = [None]          # GB ou PB, se o CSV tiver os dois
_CHECKED: set = set()
USED_MODELS: set = set()   # 'GB' e/ou 'PB' encontrados no CSV
LIGANDS: set = set()    # rótulos de resíduos do ligante (prefixo L: do gmx_MMPBSA)


def _energy(path: Path) -> pd.DataFrame:
    df = parse_energy_csv(path, MODEL[0])
    if df is None or df.empty:
        raise DataError(_("{path}: não encontrei a tabela 'Delta Energy Terms' "
                          "(rode o gmx_MMPBSA com -eo)", path=path))
    if path not in _CHECKED:
        _CHECKED.add(path)
        gap = (df["TOTAL"] - df[["VDW", "EEL", "POLAR", "NONPOLAR"]].sum(axis=1)).abs().mean()
        if gap > 0.05:
            print(_("  ⚠ {name}: ΔG total difere em {gap:.2f} kcal/mol (média) da soma vdW+elet.+polar+apolar "
                    "(termos de ligação/1-4 não nulos no delta?).", name=path.name, gap=gap))
    return df


def load_xvg(path: Path) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Energia por frame (kcal/mol) do FINAL_RESULTS_MMPBSA.csv; tempo em ns."""
    df = _energy(path)
    t_ns = df["frame"].to_numpy(float) * (NS_PER_FRAME[0] or 1.0)
    comps = {"vdw": df["VDW"].to_numpy(), "electrostatic": df["EEL"].to_numpy(),
             "polar": df["POLAR"].to_numpy(), "nonpolar": df["NONPOLAR"].to_numpy(),
             "total": df["TOTAL"].to_numpy(),
             "gas": df["GAS"].to_numpy(), "solv": df["SOLV"].to_numpy()}
    return t_ns, comps


def _res_label(resname: str, num) -> str:
    name = resname.upper()
    if name.startswith("<") or not name.strip("<>0123456789"):
        name = "LIG"                       # ligante sem nome (ex.: <0>)
    return f"{name}-{int(num)}"


def load_residues(path: Path, window: Optional[Tuple[float, float]] = None) -> pd.DataFrame:
    """Decomposição por resíduo: média e desvio (kcal/mol); aceita .csv (por frame) ou .dat (agregado).

    Com `window` (ns), o .csv usa só os frames dentro da faixa (tempo = Frame # × NS_PER_FRAME,
    como na energia vs. tempo). O .dat já vem agregado e não pode ser cortado.
    """
    cols5 = (("vdw", "VDW"), ("electrostatic", "EEL"), ("polar", "POLAR"),
             ("nonpolar", "NONPOLAR"), ("total", "TOTAL"))
    if path.suffix.lower() == ".dat":
        a = parse_decomp_dat(path)
        if a is None or a.empty:
            raise DataError(_("{path}: não encontrei 'Total Energy Decomposition' "
                              "na seção DELTAS (rode com -deo e &decomp)", path=path))
    else:
        d = parse_decomp_csv(path)
        if d is None or d.empty:
            raise DataError(_("{path}: não encontrei 'Total Energy Decomposition' "
                              "na seção DELTAS (rode com -deo e &decomp)", path=path))
        if window:
            t = d["frame"].to_numpy(float) * (NS_PER_FRAME[0] or 1.0)
            d = d[(t >= window[0]) & (t <= window[1])]
            if d.empty:
                raise DataError(_("{path}: nenhum frame da decomposição na faixa "
                                  "{a:g}–{b:g} ns", path=path, a=window[0], b=window[1]))
        g = d.groupby(["mol", "chain", "resname", "resnum"], sort=False)
        stats = g[[c for _k, c in cols5]].agg(["mean", "std"]).reset_index()
        a = pd.DataFrame({"mol": stats["mol"], "chain": stats["chain"],
                          "resname": stats["resname"], "resnum": stats["resnum"]})
        for _k, c in cols5:
            a[f"m_{c}"] = stats[(c, "mean")].to_numpy()
            a[f"sd_{c}"] = stats[(c, "std")].fillna(0).to_numpy()
    a["Residue"] = [_res_label(n, r) for n, r in zip(a["resname"], a["resnum"])]
    # mesmo rótulo em cadeias diferentes: acrescenta a cadeia
    dup = a["Residue"].duplicated(keep=False)
    if dup.any():
        a.loc[dup, "Residue"] = a.loc[dup, "Residue"] + "_" + a.loc[dup, "chain"]
    res = pd.DataFrame({"Residue": a["Residue"].to_numpy()})
    for comp, c in cols5:
        m, sd = RES_COLS[comp]
        res[m] = a[f"m_{c}"].to_numpy()
        res[sd] = a[f"sd_{c}"].to_numpy()
    LIGANDS.update(a.loc[a["mol"] == "L", "Residue"])
    return res


def norm_component(name: str) -> Optional[str]:
    n = name.lower().replace(" ", "")
    if "apolar" in n or ("non" in n and "polar" in n):
        return "nonpolar"
    if "polar" in n:
        return "polar"
    if "elec" in n:
        return "electrostatic"
    if "vdw" in n or "vanderwaals" in n:
        return "vdw"
    if "total" in n or "binding" in n:
        return "total"
    return None


def load_summary(path: Path) -> Dict[str, Tuple[float, float]]:
    """Média ± desvio entre frames dos componentes (kcal/mol)."""
    df = _energy(path)
    return {k: (float(df[c].mean()), float(df[c].std(ddof=1)))
            for k, c in (("vdw", "VDW"), ("electrostatic", "EEL"), ("polar", "POLAR"),
                         ("nonpolar", "NONPOLAR"), ("total", "TOTAL"),
                         ("gas", "GAS"), ("solv", "SOLV"))}


# ============================================================
# RESÍDUOS: LIGANTE, RANKING E TABELAS
# ============================================================

def split_label(label: str) -> Tuple[str, Optional[int]]:
    m = re.match(r"^(.*)-(\d+)$", label)
    return (m.group(1), int(m.group(2))) if m else (label, None)


def detect_ligand(df: pd.DataFrame) -> Optional[str]:
    known = [r for r in df["Residue"] if r in LIGANDS]
    if known:
        return known[-1]
    cands = [r for r in df["Residue"]
             if split_label(r)[0].upper() not in STANDARD_RESIDUES]
    return cands[-1] if cands else None


def make_ranking(df: pd.DataFrame, drop: Optional[str]) -> pd.DataFrame:
    d = df[df["Residue"] != drop].copy() if drop else df.copy()
    d = d.sort_values("total", ascending=True).reset_index(drop=True)
    d.insert(0, "Rank", range(1, len(d) + 1))
    return d


def pm_table(rank: pd.DataFrame, label: str, lang: str) -> pd.DataFrame:
    """Tabela 'média ± desvio' (cabeçalhos no idioma escolhido)."""
    S = REPORT_TEXT[lang]
    t = pd.DataFrame({S["rank"]: rank["Rank"], S["complex"]: label,
                      S["residue"]: rank["Residue"]})
    for comp in COMPONENTS:
        m, s = RES_COLS[comp]
        t[S["pm_cols"][comp]] = [f"{a:.2f} ± {b:.2f}"
                                 for a, b in zip(rank[m], rank[s])]
    return t


def rank_text(rank: pd.DataFrame, cfg: configparser.ConfigParser,
              lang: str) -> str:
    S = REPORT_TEXT[lang]
    cols, heads = ["Rank", "Residue"], [S["rank"], S["residue"]]
    for comp in COMPONENTS:
        cols += list(RES_COLS[comp])
        heads += [cfg[f"text_{lang}"][f"comp_{comp}"], S["sd"]]
    sub = rank[cols].copy()
    sub.columns = heads
    return sub.to_string(index=False, float_format=lambda x: f"{x:.2f}")


# ============================================================
# SAÍDA DE FIGURAS
# ============================================================

def finish(ctx: Ctx, fig, stem: str, group: str, rect=None) -> None:
    """Mostra (opcional), salva nos formatos escolhidos e fecha a figura."""
    if rect:
        fig.tight_layout(rect=rect)
    else:
        fig.tight_layout()
    if ctx.interactive:
        if matplotlib.get_backend().lower() == "agg":
            ctx.report.warn(_("Sem backend gráfico disponível; "
                              "janela interativa ignorada."))
        else:
            print(_("  ... ajuste '{stem}' na janela e feche-a para salvar.", stem=stem))
            plt.show()
    dpi = cfg_float(ctx.cfg, "style", "dpi")
    for fmt in ctx.formats:
        path = ctx.target_dir(group, fmt) / f"{pfx()}{stem}.{fmt}"
        if fmt == "png":
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
        else:
            fig.savefig(path, bbox_inches="tight")
        ctx.report.file(path)
    plt.close(fig)


DECIMAL_COMMA = [False]     # True quando o idioma é português


class CommaFormatter(ScalarFormatter):
    """ScalarFormatter com vírgula decimal (padrão brasileiro)."""

    def __call__(self, x, pos=None):
        return super().__call__(x, pos).replace(".", ",")


def style_axes(ax, minor_x: bool = True) -> None:
    """Acabamento comum: moldura, ticks para dentro com ticks menores."""
    for sp in ax.spines.values():
        sp.set_linewidth(1.2)
    ax.set_axisbelow(True)
    if DECIMAL_COMMA[0]:
        ax.yaxis.set_major_formatter(CommaFormatter())
        if minor_x:                      # eixo x numérico (não categórico)
            ax.xaxis.set_major_formatter(CommaFormatter())
    ax.tick_params(axis="both", which="major", direction="in", length=5,
                   width=1.1, pad=5, top=False, right=False)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis="y", which="minor", direction="in", length=2.5,
                   width=0.8, right=False)
    if minor_x:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(axis="x", which="minor", direction="in", length=2.5,
                       width=0.8, top=False)


def legend_on_top(ax, n_items: int, ncol: Optional[int] = None) -> None:
    """Legenda em linha no topo, numa faixa livre acima dos dados."""
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.16 * (hi - lo))
    ax.legend(frameon=False, loc="upper center", ncol=ncol or min(n_items, 3),
              columnspacing=1.8, handlelength=1.8, borderaxespad=0.4)


# ============================================================
# GRÁFICOS
# ============================================================

def plot_time(ctx: Ctx, key: str, systems: List[System],
              data: Dict[str, Tuple[np.ndarray, Dict[str, np.ndarray]]],
              window: Tuple[float, float]) -> None:
    lw = cfg_float(ctx.cfg, "style", "line_width")
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    drawn = 0
    for s in systems:
        if s.key not in data:
            continue
        t, comps = data[s.key]
        m = (t >= window[0]) & (t <= window[1])
        if not m.any():
            ctx.report.warn(_("{label}: nenhum frame entre {a} e {b} ns "
                              "(dados: {lo:.1f}–{hi:.1f} ns).", label=s.label,
                              a=window[0], b=window[1], lo=t.min(), hi=t.max()))
            continue
        ax.plot(t[m], comps[key][m], label=s.label, color=s.color,
                linewidth=lw, alpha=0.9)
        drawn += 1
    if not drawn:
        plt.close(fig)
        ctx.report.warn(_("Gráfico '{key}' não gerado: sem dados na faixa.", key=key))
        return
    ax.set_xlabel(ctx.t("time_xlabel"), fontweight="bold", labelpad=8)
    ax.set_ylabel(ctx.t("energy_ylabel"), fontweight="bold", labelpad=8)
    ax.set_title(ctx.t(f"title_{key}"), fontweight="bold", pad=14)
    ax.set_xlim(*window)
    ax.axhline(0, color="black", linewidth=0.7, linestyle="--", alpha=0.6)
    ax.grid(False)
    style_axes(ax)
    legend_on_top(ax, drawn)
    finish(ctx, fig, f"energy_{key}_vs_time", "time")


def plot_summary(ctx: Ctx, systems: List[System],
                 sums: Dict[str, Dict[str, Tuple[float, float]]]) -> None:
    sel = [s for s in systems if s.key in sums]
    if not sel:
        return
    n = len(sel)
    width = 0.8 / n
    x = np.arange(len(ECOMPS))
    fig, ax = plt.subplots(figsize=(max(7.2, 1.6 * len(ECOMPS) + n), 5.2))
    for i, s in enumerate(sel):
        means = [sums[s.key].get(c, (np.nan, np.nan))[0] for c in ECOMPS]
        sds = [sums[s.key].get(c, (np.nan, np.nan))[1] for c in ECOMPS]
        pos = x + (i - (n - 1) / 2) * width
        ax.bar(pos, means, width=width, color=s.color, label=s.label,
               edgecolor="none", zorder=3)
        ax.errorbar(pos, means, yerr=sds, fmt="none", ecolor="#444444",
                    elinewidth=1.1, capsize=3, capthick=1.1, zorder=4)
    ax.axhline(0, color="black", linewidth=0.8, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([ctx.t(f"comp_{c}") for c in ECOMPS],
                       fontweight="bold")
    ax.set_ylabel(ctx.t("energy_ylabel"), fontweight="bold", labelpad=8)
    ax.set_title(ctx.t("summary_title"), fontweight="bold", pad=14)
    ax.grid(axis="y", linewidth=0.6, alpha=0.3, zorder=0)
    style_axes(ax, minor_x=False)
    legend_on_top(ax, n)
    finish(ctx, fig, "energy_components_mean", "summary")


def plot_residues_top(ctx: Ctx, s: System, rank: pd.DataFrame,
                      top_n: int, comp_colors: List[str]) -> None:
    top = rank.head(top_n)
    residues = top["Residue"].tolist()
    k = len(residues)
    x = np.arange(k)
    width = 0.15
    fig, ax = plt.subplots(figsize=(min(11.0, max(7.5, 0.95 * k + 3.0)), 5.4))
    for i, comp in enumerate(COMPONENTS):
        m, sd = RES_COLS[comp]
        pos = x + (i - 2) * width
        ax.bar(pos, top[m].to_numpy(), width=width, color=comp_colors[i],
               label=ctx.t(f"comp_{comp}"), edgecolor="none", zorder=3)
        ax.errorbar(pos, top[m].to_numpy(), yerr=top[sd].to_numpy(),
                    fmt="none", ecolor="#444444", elinewidth=1.1,
                    capsize=3, capthick=1.1, zorder=4)
    ax.axhline(0, color="#c8a800", linewidth=1.2, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(residues, rotation=45, ha="right", fontweight="bold")
    ax.set_xlim(-0.6, k - 0.4)
    ax.set_xlabel(ctx.t("residue_xlabel"), fontweight="bold", labelpad=8)
    ax.set_ylabel(ctx.t("energy_ylabel"), fontweight="bold", labelpad=8)
    # título em negrito + subtítulo (nome do sistema) em fonte normal
    ax.set_title(ctx.t("residue_title"), fontweight="bold",
                 fontsize=cfg_float(ctx.cfg, "font", "title_size") + 2, pad=30)
    ax.annotate(s.label, xy=(0.5, 1.0), xycoords="axes fraction",
                xytext=(0, 8), textcoords="offset points", ha="center",
                va="bottom", fontsize=cfg_float(ctx.cfg, "font", "label_size"))
    ax.grid(axis="y", linewidth=0.6, alpha=0.3, zorder=0)
    style_axes(ax, minor_x=False)
    legend_on_top(ax, len(COMPONENTS), ncol=len(COMPONENTS))
    finish(ctx, fig, f"{s.slug}_per_residue_top{top_n}", "residues")


def plot_residue_profile(ctx: Ctx, s: System, rank: pd.DataFrame) -> None:
    d = rank.copy()
    nums = [split_label(r)[1] for r in d["Residue"]]
    if any(n is None for n in nums):
        nums = list(range(1, len(d) + 1))
    d["_n"] = nums
    d = d.sort_values("_n")
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.bar(d["_n"], d["total"], width=1.0, color=s.color, edgecolor="none",
           zorder=3)
    ax.axhline(0, color="black", linewidth=0.8, zorder=2)
    for i, (_i, row) in enumerate(rank.head(5).iterrows()):
        n = split_label(row["Residue"])[1]
        if n is None:
            continue
        # alturas escalonadas + linha de ligação, para vizinhos não se sobreporem
        ax.annotate(row["Residue"], (n, row["total"]),
                    xytext=(0, -14 - 14 * i), textcoords="offset points",
                    ha="center", va="top", fontsize=9,
                    arrowprops=dict(arrowstyle="-", linewidth=0.6,
                                    color="#555555"))
    ax.set_xlabel(ctx.t("residue_number_xlabel"), fontweight="bold",
                  labelpad=8)
    ax.set_ylabel(ctx.t("energy_ylabel"), fontweight="bold", labelpad=8)
    ax.set_title(f"{ctx.t('profile_title')} — {s.label}", fontweight="bold",
                 pad=14)
    ax.margins(x=0.01)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo - 0.50 * (hi - lo), hi)    # espaço para os rótulos dos melhores
    ax.grid(axis="y", linewidth=0.6, alpha=0.3, zorder=0)
    style_axes(ax)
    finish(ctx, fig, f"{s.slug}_per_residue_profile", "residues")


# ============================================================
# MENUS
# ============================================================

PLOT_ITEMS = [
    ("total", "Energia no tempo — Total"),
    ("vdw", "Energia no tempo — van der Waals"),
    ("electrostatic", "Energia no tempo — Eletrostática"),
    ("polar", "Energia no tempo — Solvatação polar"),
    ("nonpolar", "Energia no tempo — Solvatação apolar"),
    ("summary", "Médias ± desvio dos componentes (FINAL_RESULTS_MMPBSA.csv)"),
    ("residues", "Decomposição por resíduo (FINAL_DECOMP_MMPBSA.csv)"),
]
TIME_KEYS = ["total", "vdw", "electrostatic", "polar", "nonpolar"]
ECOMPS = list(COMPONENTS)      # componentes de energia (barras de média e relatório)


def enable_gbsa() -> None:
    """MM/GBSA: acrescenta ΔG_gas (GGAS) e ΔG_solv (GSOLV), como nos artigos de GBSA."""
    if "gas" in ECOMPS:
        return
    ECOMPS[:] = ["vdw", "electrostatic", "polar", "nonpolar", "gas", "solv", "total"]
    TIME_KEYS.extend(["gas", "solv"])
    k = [i for i, (key, _l) in enumerate(PLOT_ITEMS) if key == "nonpolar"][0] + 1
    PLOT_ITEMS[k:k] = [("gas", "Energia no tempo — Fase gasosa (GGAS)"),
                       ("solv", "Energia no tempo — Solvatação (GSOLV)")]


def available_plots(systems: List[System]) -> List[Tuple[str, str]]:
    has_xvg = any(s.xvg for s in systems)
    has_sum = any(s.summary for s in systems)
    has_res = any(s.res for s in systems)
    out = []
    for key, label in PLOT_ITEMS:
        if key in TIME_KEYS and has_xvg:
            out.append((key, label))
        elif key == "summary" and has_sum:
            out.append((key, label))
        elif key == "residues" and has_res:
            out.append((key, label))
    return out


def choose_plots(items: List[Tuple[str, str]], preset: Optional[str]) -> List[str]:
    keys = [k for k, _k in items]

    def parse(text: str) -> List[str]:
        t = text.strip().lower()
        if t in ("all", "todos", "a", "*"):
            return keys
        parts = [p for p in re.split(r"[\s,;]+", t) if p]
        if parts and all(p in keys for p in parts):     # por nome (--plots)
            return [k for k in keys if k in parts]
        return [keys[i] for i in parse_indices(text, len(keys))]

    if preset is None and not AUTO:
        for i, (_k, label) in enumerate(items, 1):
            print(f"  [{i}]  {_(label)}")
        print(_("  [all]  todos") + "\n")
    return ask("Quais gráficos? (ex.: 1,3 ou all)", parse, preset)


def fmt_systems(systems: List[System]) -> None:
    if AUTO:
        return
    width = max(len(s.key) for s in systems) + 2
    print(_("  {n} sistema(s) encontrado(s):", n=len(systems)) + "\n")
    for i, s in enumerate(systems, 1):
        have = [_(n) for n, f in (("energia vs. tempo", s.xvg),
                                  ("resíduos", s.res),
                                  ("médias", s.summary)) if f]
        print(f"  [{i}]  {s.key:<{width}}{dim(' · '.join(have))}")
    print()


# ============================================================
# EXECUÇÃO
# ============================================================

def run(args, report: Report) -> None:
    global AUTO
    AUTO = bool(args.yes)
    MODEL[0] = args.model
    root = Path(args.root).resolve()
    outdir = Path(args.out).resolve() if args.out else root / OUT_DIR_NAME
    cfg_path = Path(args.config).resolve() if args.config else root / CONFIG_NAME

    # ---------- idioma: interface, gráficos, tabelas e relatório ----------
    print()
    print(bold("  gmx_mmpbsa_plot"))
    if args.lang is None and not AUTO:
        print(dim("  Idioma / Language") + "\n")
        print("  [1]  Português (Brasil)")
        print("  [2]  English (US)\n")
    lang = ask("Idioma / Language", parse_lang, args.lang,
               default=peek_saved_lang(cfg_path))
    UI_LANG[0] = lang
    DECIMAL_COMMA[0] = lang == "pt"

    if not root.is_dir():
        raise UsageError(_("Pasta não encontrada: {root}", root=root))
    print(dim("  " + _("gráficos de MM/PBSA e MM/GBSA (gmx_MMPBSA / Amber)")))
    print(dim(_("  Pasta raiz: {root}", root=root)))
    report.base = outdir

    systems = discover(root, outdir)
    if not systems:
        raise UsageError(_("Nenhum arquivo ({a}, {b}, {c}) encontrado em {root} "
                           "ou nas subpastas.", a=XVG_NAME, b=RES_NAME,
                           c=SUM_NAME, root=root))

    cfg, cfg_changed = load_config(cfg_path)
    apply_style(cfg)
    fam = cfg["font"].get("family", "Arial")
    font_used = resolve_font(fam)
    if font_used != fam:
        print(dim(_("  Fonte '{fam}' não instalada; usando '{used}'.",
                    fam=fam, used=font_used)))

    # ---------- 1. sistemas ----------
    step("Sistemas")
    fmt_systems(systems)
    idx = ask("Quais sistemas? (ex.: 1,3 ou all)",
              lambda t: parse_indices(t, len(systems)), args.systems,
              default="all")
    chosen = [systems[i] for i in idx]

    # ---------- 2. legendas (padrão = nome da pasta; editável) ----------
    if not AUTO:
        step("Legendas")
        print(dim(_("  Enter mantém o nome entre colchetes. Os nomes ficam em "
                    "{cfg}\n  e podem ser editados lá depois.", cfg=cfg_path.name))
              + "\n")
    for s in chosen:
        current = cfg["names"].get(s.key, s.key)
        name = ask(s.key, parse_name, None, default=current)
        s.label = name
        if cfg["names"].get(s.key) != name:
            cfg["names"][s.key] = name
            cfg_changed = True

    # ---------- modelo de solvatação: GBSA ou PBSA ----------
    en_files = [s.xvg for s in chosen if s.xvg]
    if en_files:
        present = set().union(*(detect_models(f) for f in en_files))
        if present:
            step("Modelo")
            if args.model is None:
                print("  [1]  MM/GBSA   " + dim("(Generalized Born)"
                      + ("" if "GB" in present else "  — " + _("não encontrado nos arquivos"))))
                print("  [2]  MM/PBSA   " + dim("(Poisson-Boltzmann)"
                      + ("" if "PB" in present else "  — " + _("não encontrado nos arquivos"))) + "\n")
            default = {"GB": "1", "PB": "2"}[sorted(present)[0]]
            MODEL[0] = ask("Qual método foi usado? (1=MM/GBSA, 2=MM/PBSA)", parse_model, args.model, default=default)
            if MODEL[0] not in present:
                name = {"GB": "MM/GBSA", "PB": "MM/PBSA"}[MODEL[0]]
                raise UsageError(_("Os arquivos não têm resultados {name}; "
                                   "encontrado: {found}.", name=name,
                                   found=", ".join({"GB": "MM/GBSA", "PB": "MM/PBSA"}[m] for m in sorted(present))))
            if not AUTO:
                print()
    if MODEL[0] == "GB":
        enable_gbsa()

    # ---------- 3. gráficos ----------
    items = available_plots(chosen)
    if not items:
        raise UsageError(_("Os sistemas escolhidos não têm arquivos plotáveis."))
    step("Gráficos")
    plots = choose_plots(items, args.plots)
    do_time = [k for k in TIME_KEYS if k in plots]
    do_summary = "summary" in plots
    do_res = "residues" in plots

    # ---------- 4. parâmetros específicos ----------
    window: Optional[Tuple[float, float]] = None
    if do_time or do_res:
        step("Parâmetros")
    if do_time or do_res:
        NS_PER_FRAME[0] = ask("Tempo em ns de cada frame da trajetória (coluna 'Frame #' do CSV)",
                              parse_pos_float, args.ns_per_frame)
        spans = []
        for s in chosen:
            if s.xvg:
                try:
                    t, _unused = load_xvg(s.xvg)
                    spans.append((t.min(), t.max()))
                except DataError:
                    pass
        hint = (_(" (dados: {lo:g}–{hi:g} ns)",
                  lo=min(a for a, _b in spans),
                  hi=max(b for _a, b in spans))) if spans else ""
        if args.time:
            preset = f"{args.time[0]} {args.time[1]}"
        else:
            preset = None
        window = ask(_("Faixa de tempo em ns, início e fim{hint}", hint=hint),
                     parse_window, preset)
        if not AUTO:
            print()

    top_n: Optional[int] = 5
    keep_ligand = False
    comp_colors = list(COMPONENT_COLORS)
    if do_res:
        if args.top is None and not AUTO:
            print(_("  Decomposição por resíduo — quantos mostrar?") + "\n")
            print("  [1]  Top 5")
            print("  [2]  Top 10")
            print("  [3]  " + _("Todos") + "   "
                  + dim(_("(não recomendado: com centenas de resíduos o gráfico de\n"
                          "                barras fica ilegível; gera só o perfil por posição\n"
                          "                na sequência e o CSV completo)")) + "\n")
        top_n = ask("Quantos resíduos?", parse_top, args.top, default=None)
        if top_n is None and not AUTO and args.top is None:
            print()
            if not ask("'Todos' não é recomendado. Continuar mesmo assim?",
                       parse_yn, None, default="n"):
                top_n = ask("Quantos resíduos? (1=Top 5, 2=Top 10)",
                            parse_top_limited, None)
        if not AUTO:
            print()

        # detecção do ligante
        det = {}
        for s in chosen:
            if s.res:
                try:
                    s.ligand = detect_ligand(load_residues(s.res))
                except DataError:
                    pass
                det[s.key] = s.ligand
        if not AUTO:
            print(_("  Ligante detectado em cada sistema:") + "\n")
            for s in chosen:
                if s.res:
                    print(f"    {s.label:<28}{s.ligand or dim(_('não detectado'))}")
            print()
        if any(det.values()):
            if args.ligand:
                keep_ligand = args.ligand == "keep"
            else:
                keep_ligand = ask(
                    "Manter o ligante no ranking? (mantido, ele tende a "
                    "dominar o topo)", parse_yn, None, default="n")
                print()
        for s in chosen:
            s.drop_ligand = None if keep_ligand else s.ligand

    # ---------- 5. cores ----------
    user_colors: List[str] = []
    if (do_time or do_summary or do_res) and not AUTO:
        step("Cores")
        print(dim(_("  Hex separados por vírgula (ex.: #0072B2,#D55E00). "
                    "Enter usa a paleta padrão.")) + "\n")
    if args.colors:
        user_colors = ask("Cores dos sistemas, na ordem escolhida",
                          parse_hex_list, args.colors)
    elif not AUTO and (do_time or do_summary):
        user_colors = ask("Cores dos sistemas, na ordem escolhida",
                          lambda t: parse_hex_list(t) if t.strip() else [],
                          None, default="")
    if user_colors and len(user_colors) < len(chosen):
        report.warn(_("Menos cores do que sistemas; o restante usa a paleta padrão."))
    for i, s in enumerate(chosen):
        s.color = user_colors[i] if i < len(user_colors) else PALETTE[i % len(PALETTE)]

    if do_res:
        if args.comp_colors:
            comp_colors = ask("Cores dos 5 componentes (vDW, elet., polar, apolar, total)",
                              parse_hex_list, args.comp_colors)
        elif not AUTO:
            txt = ask("Cores dos 5 componentes (vDW, elet., polar, apolar, total)",
                      lambda t: parse_hex_list(t) if t.strip() else [], None,
                      default="")
            if txt:
                comp_colors = txt
        if len(comp_colors) < 5:
            report.warn(_("Menos de 5 cores de componentes; completando com o padrão."))
            comp_colors += COMPONENT_COLORS[len(comp_colors):]

    # ---------- 6. saída ----------
    step("Saída")
    if cfg["general"].get("language") != lang:
        cfg["general"]["language"] = lang
        cfg_changed = True
    if args.layout is None and not AUTO:
        names = GROUP_DIRS[lang]
        print(_("  Organização dos arquivos gerados:") + "\n")
        print("  [1]  " + _("por tipo de gráfico (recomendado)").ljust(34)
              + dim(f"{names['time']}/, {names['residues']}/, {names['summary']}/"))
        print("  [2]  " + _("por formato").ljust(34) + dim("png/, pdf/, svg/"))
        print("  [3]  " + _("tudo na mesma pasta") + "\n")
    layout = ask("Organização", parse_layout, args.layout, default="1")
    if not AUTO:
        print()
    formats = ask("Formatos (png, pdf, svg; Enter = todos)",
                  parse_formats, args.formats, default="all")
    export_csv = False
    if do_time:
        export_csv = ask("Exportar também os dados de energia em CSV "
                         "(ns e kcal/mol)?", parse_yn,
                         "s" if args.csv else None, default="n")
    interactive = ask("Abrir janela interativa para ajustar cada figura antes "
                      "de salvar?", parse_yn,
                      "s" if args.interactive else None, default="n")

    outdir.mkdir(parents=True, exist_ok=True)
    if cfg_changed:
        save_config(cfg, cfg_path)

    ctx = Ctx(root=root, outdir=outdir, cfg=cfg, report=report, formats=formats,
              interactive=interactive, export_csv=export_csv,
              lang=lang, layout=layout)
    S = REPORT_TEXT[lang]

    # ---------- relatório: cabeçalho e escolhas ----------
    report.add(S["title"])
    report.add(f"{S['date']}: {dt.datetime.now():%Y-%m-%d %H:%M:%S}")
    report.add(f"{S['root']}: {root}")
    report.add(f"{S['out']}: {outdir}")
    report.add(f"{S['cfg']}: {cfg_path}")
    report.section(S["settings"])
    report.add(f"{S['lang']}: {'English' if lang == 'en' else 'Português'}")
    report.add(f"{S['font']}: {font_used}"
               + ("" if font_used == cfg["font"].get("family", "Arial")
                  else f" ({S['font_sub']} {cfg['font'].get('family', 'Arial')})"))
    report.add(f"{S['systems']}:")
    for s in chosen:
        report.add(f"  {s.key:<32} -> {s.label}   ({S['color']} {s.color})")
    if MODEL[0]:
        report.add(f"Modelo / Model: {'MM/GBSA (Generalized Born)' if MODEL[0] == 'GB' else 'MM/PBSA (Poisson-Boltzmann)'}")
    report.add(f"{S['plots']}: {', '.join(plots)}")
    if window:
        report.add(f"{S['window']}: {window[0]:g} – {window[1]:g} ns")
    if do_res:
        top_txt = S["all"] if top_n is None else f"Top {top_n}"
        report.add(f"{S['residues']}: {top_txt}")
        report.add(f"{S['ligand_kept']}: {S['yes'] if keep_ligand else S['no']}")
        for s in chosen:
            if s.res:
                report.add(f"  {S['ligand_det']} {s.label}: "
                           f"{s.ligand or S['none']}")
    report.add(f"{S['formats']}: {', '.join(formats)}")
    report.add(f"{S['layout']}: {S['layout_names'][layout]}")
    report.add(S["units"])

    print()
    print(bold(_("Gerando")) + dim(f"  →  {outdir}"))
    print(dim("─" * 64))

    # ---------- energia no tempo ----------
    if do_time:
        print(bold("\n  " + _("Energia no tempo")))
        data = {}
        for s in chosen:
            if not s.xvg:
                report.warn(_("{label}: sem {fname}; ignorado nos gráficos de tempo.",
                              label=s.label, fname=XVG_NAME))
                continue
            try:
                data[s.key] = load_xvg(s.xvg)
            except DataError as e:
                report.warn(str(e))
        for key in do_time:
            plot_time(ctx, key, chosen, data, window)

        report.section(S["time_sec"].format(a=window[0], b=window[1]))
        hdr = f"{S['system']:<28}{S['frames']:>7}" + "".join(
            f"{ctx.t('comp_' + c):>20}" for c in ECOMPS)
        report.add(hdr)
        for s in chosen:
            if s.key not in data:
                continue
            t, comps = data[s.key]
            m = (t >= window[0]) & (t <= window[1])
            if not m.any():
                continue
            row = f"{s.label[:27]:<28}{int(m.sum()):>7}"
            for c in ECOMPS:
                v = comps[c][m]
                sd = v.std(ddof=1) if v.size > 1 else 0.0
                row += f"{v.mean():>12.2f} ± {sd:<5.2f}"
            report.add(row)
            if ctx.export_csv:
                df = pd.DataFrame({"time_ns": t[m],
                                   **{c: comps[c][m] for c in ECOMPS}})
                p = ctx.target_dir("tables") / f"{pfx()}{s.slug}_energy_vs_time.csv"
                df.to_csv(p, index=False)
                report.file(p)

    # ---------- médias do summary ----------
    if do_summary:
        print(bold("\n  " + _("Médias dos componentes")))
        sums = {}
        for s in chosen:
            if not s.summary:
                report.warn(_("{label}: sem {fname}; ignorado nas médias.",
                              label=s.label, fname=SUM_NAME))
                continue
            try:
                sums[s.key] = load_summary(s.summary)
            except DataError as e:
                report.warn(str(e))
        plot_summary(ctx, chosen, sums)
        report.section(S["avg_sec"])
        report.add(f"{S['system']:<28}" + "".join(
            f"{ctx.t('comp_' + c):>20}" for c in ECOMPS))
        for s in chosen:
            if s.key in sums:
                row = f"{s.label[:27]:<28}"
                for c in ECOMPS:
                    m, sd = sums[s.key].get(c, (float('nan'), float('nan')))
                    row += f"{m:>12.2f} ± {sd:<5.2f}"
                report.add(row)

    # ---------- resíduos ----------
    if do_res:
        print(bold("\n  " + _("Decomposição por resíduo")))
        report.section(S["res_sec"])
        for s in chosen:
            if not s.res:
                report.warn(_("{label}: sem {fname}; ignorado nos resíduos.",
                              label=s.label, fname=RES_NAME))
                continue
            try:
                use_win = window if s.res.suffix.lower() != ".dat" else None
                df = load_residues(s.res, use_win)
                if use_win:
                    report.add(f"[{s.label}] " + S["res_window"].format(
                        a=window[0], b=window[1]))
                elif window:
                    report.warn(_("{label}: {fname} já vem agregado (média ± DP); resíduos "
                                  "calculados na trajetória inteira, não na faixa "
                                  "{a:g}–{b:g} ns.", label=s.label, fname=s.res.name,
                                  a=window[0], b=window[1]))
            except DataError as e:
                report.warn(str(e))
                continue
            print(f"  {s.label}")
            rank = make_ranking(df, s.drop_ligand)
            full = ctx.target_dir("tables") / f"{pfx()}{s.slug}_residue_ranking_complete.csv"
            rank.to_csv(full, index=False)
            report.file(full)
            if top_n is None:
                plot_residue_profile(ctx, s, rank)
                n_show = 10
                head = S["full_saved"]
            else:
                table = pm_table(rank.head(top_n), s.label, lang)
                tp = ctx.target_dir("tables") / f"{pfx()}{s.slug}_top_{top_n}_table.csv"
                table.to_csv(tp, index=False)
                report.file(tp)
                plot_residues_top(ctx, s, rank, top_n, comp_colors)
                n_show = top_n
                head = f"Top {top_n}"
            removed = (" — " + S["ligand_excl"].format(l=s.drop_ligand)
                       if s.drop_ligand else " — " + S["ligand_keep"])
            report.add()
            report.add(f"[{s.label}] {head}{removed}")
            report.add(rank_text(rank.head(n_show), cfg, lang))

    # ---------- rodapé ----------
    report.section(S["files"])
    for f in report.files:
        report.add(str(f))
    if report.warnings:
        report.section(S["warnings"])
        for w in report.warnings:
            report.add(f"- {w}")


def detect_models(path: Path) -> set:
    """Modelos de solvatação presentes no CSV: {'GB'}, {'PB'} ou os dois."""
    found = set()
    for ln in path.read_text(errors="replace").splitlines():
        low = ln.lower()
        if "generalized born" in low:
            found.add("GB")
        elif "poisson" in low:
            found.add("PB")
    return found


def parse_model(text: str) -> str:
    t = text.strip().lower().replace("/", "").replace(" ", "")
    if t in ("1", "gb", "gbsa", "mmgbsa", "gbm"):
        return "GB"
    if t in ("2", "pb", "pbsa", "mmpbsa"):
        return "PB"
    raise ValueError(_("escolha 1 (MM/GBSA) ou 2 (MM/PBSA)"))


def pfx() -> str:
    """Prefixo dos arquivos conforme o método: GB_ ou PB_ (vazio se não souber)."""
    return f"{MODEL[0]}_" if MODEL[0] else ""


def method_name() -> str:
    """MM/GBSA, MM/PBSA ou os dois, conforme o modelo lido do CSV."""
    if USED_MODELS == {"GB"}:
        return "MM/GBSA"
    if USED_MODELS == {"PB"}:
        return "MM/PBSA"
    return "MM/PBSA" if not USED_MODELS else "MM/PBSA + MM/GBSA"


def write_report(report: Report, outdir: Path, status: str) -> Optional[Path]:
    if not report.lines:
        return None
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{pfx()}mmpbsa_report_{dt.datetime.now():%Y%m%d_%H%M%S}.txt"
    lines = [ln.replace("MM/PBSA", method_name()) for ln in report.lines]
    if status != "ok":
        lines.insert(1, f"STATUS: {status}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="MM/PBSA and MM/GBSA plots for gmx_MMPBSA (Amber) output, with an "
                    "interactive menu. Any option below skips the "
                    "corresponding menu question.")
    p.add_argument("root", nargs="?", default=".",
                   help="root folder containing one subfolder per system "
                        "(default: .)")
    p.add_argument("-o", "--out",
                   help=f"output folder (default: a new <root>/{OUT_DIR_NAME}_YYYYMMDD_HHMMSS "
                        "each run, so nothing is overwritten)")
    p.add_argument("--config",
                   help=f"configuration file (default: <root>/{CONFIG_NAME})")
    p.add_argument("--lang", metavar="{pt,en}",
                   help="language of the interface, plots, tables and report "
                        "(skips the initial language menu)")
    p.add_argument("--systems", help="system indices, e.g. 1,3 or all")
    p.add_argument("--plots",
                   help="plots by name or number: total, vdw, electrostatic, "
                        "polar, nonpolar, summary, residues, or all")
    p.add_argument("--ns-per-frame", metavar="NS",
                   help="time in ns of one trajectory frame ('Frame #' column of the CSV; "
                        "e.g. 0.01 for 10 ps)")
    p.add_argument("--model", choices=["GB", "PB"],
                   help="solvation model, if the CSV has both")
    p.add_argument("--time", nargs=2, type=float, metavar=("START", "END"),
                   help="time range in ns")
    p.add_argument("--top",
                   help="5, 10 or all (per-residue decomposition)")
    p.add_argument("--ligand", choices=["keep", "drop"],
                   help="keep or drop the ligand in the residue ranking")
    p.add_argument("--formats", help="output formats: png,pdf,svg")
    p.add_argument("--layout", metavar="{type,format,flat}",
                   help="output folder layout: by plot type, by file format, "
                        "or a single folder")
    p.add_argument("--colors",
                   help="comma-separated hex colors for the systems")
    p.add_argument("--comp-colors",
                   help="comma-separated hex colors for the 5 components "
                        "(per-residue plot)")
    p.add_argument("--csv", action="store_true",
                   help="also export the energy-vs-time data as CSV")
    p.add_argument("--interactive", action="store_true",
                   help="open an interactive window before saving each figure")
    p.add_argument("-y", "--yes", action="store_true",
                   help="never prompt; use defaults (all systems, all formats, "
                        "ligand dropped). --plots is required, plus --time "
                        "and --top when those plots are selected")
    return p


def main() -> int:
    args = build_parser().parse_args()
    report = Report()
    if args.out is None:        # padrão: pasta nova a cada execução (não sobrescreve)
        stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        args.out = str(Path(args.root).resolve() / f"{OUT_DIR_NAME}_{stamp}")
    outdir = Path(args.out).resolve()
    status = "ok"
    code = 0
    try:
        run(args, report)
    except (KeyboardInterrupt, EOFError):
        print("\n\n" + _("Interrompido."))
        status = _("INTERROMPIDO pelo usuário")
        code = 130
    except (UsageError, DataError) as e:
        print("\n" + _("Erro: {err}", err=e), file=sys.stderr)
        status = _("ERRO: {err}", err=e)
        code = 2
    except Exception:
        traceback.print_exc()
        status = _("ERRO inesperado (veja o terminal)")
        code = 1
    path = write_report(report, outdir, status)
    print()
    print(dim("─" * 64))
    if code == 0:
        n_fig = sum(1 for f in report.files if f.suffix != ".csv")
        n_tab = len(report.files) - n_fig
        print(bold(_("Concluído.")) + "  "
              + _("{f} arquivo(s) de figura, {t} tabela(s)", f=n_fig, t=n_tab)
              + (", " + _("{w} aviso(s)", w=len(report.warnings))
                 if report.warnings else "") + ".")
        print("  " + _("Saída:").ljust(11) + str(outdir))
    if path:
        print("  " + _("Relatório:").ljust(11) + path.name)
    print()
    return code


# ============================================================
# TRADUÇÃO DA INTERFACE (chave = texto em português)
# ============================================================

EN_UI: Dict[str, str] = {
    # etapas do menu
    "Sistemas": "Systems",
    "Legendas": "Labels",
    "Gráficos": "Plots",
    "Parâmetros": "Parameters",
    "Cores": "Colors",
    "Saída": "Output",
    # sistemas e gráficos
    "  {n} sistema(s) encontrado(s):": "  {n} system(s) found:",
    "energia vs. tempo": "energy vs. time",
    "resíduos": "residues",
    "médias": "averages",
    "Quais sistemas? (ex.: 1,3 ou all)": "Which systems? (e.g. 1,3 or all)",
    "  Enter mantém o nome entre colchetes. Os nomes ficam em {cfg}\n"
    "  e podem ser editados lá depois.":
        "  Enter keeps the name in brackets. Names are saved in {cfg}\n"
        "  and can be edited there later.",
    "Energia no tempo — Total": "Energy vs. time — Total",
    "Energia no tempo — van der Waals": "Energy vs. time — van der Waals",
    "Energia no tempo — Eletrostática": "Energy vs. time — Electrostatic",
    "Energia no tempo — Solvatação polar": "Energy vs. time — Polar solvation",
    "Energia no tempo — Solvatação apolar": "Energy vs. time — Nonpolar solvation",
    "Médias ± desvio dos componentes (FINAL_RESULTS_MMPBSA.csv)":
        "Component averages ± SD (FINAL_RESULTS_MMPBSA.csv)",
    "Decomposição por resíduo (FINAL_DECOMP_MMPBSA.csv)":
        "Per-residue decomposition (FINAL_DECOMP_MMPBSA.csv)",
    "  [all]  todos": "  [all]  all",
    "Quais gráficos? (ex.: 1,3 ou all)": "Which plots? (e.g. 1,3 or all)",
    "Os sistemas escolhidos não têm arquivos plotáveis.":
        "The chosen systems have no plottable files.",
    # parâmetros
    " (dados: {lo:g}–{hi:g} ns)": " (data: {lo:g}–{hi:g} ns)",
    "Faixa de tempo em ns, início e fim{hint}":
        "Time range in ns, start and end{hint}",
    "  Decomposição por resíduo — quantos mostrar?":
        "  Per-residue decomposition — how many to show?",
    "Todos": "All",
    "(não recomendado: com centenas de resíduos o gráfico de\n"
    "                barras fica ilegível; gera só o perfil por posição\n"
    "                na sequência e o CSV completo)":
        "(not recommended: with hundreds of residues the bar chart is\n"
        "                unreadable; only the sequence-position profile and\n"
        "                the complete CSV are generated)",
    "Quantos resíduos?": "How many residues?",
    "'Todos' não é recomendado. Continuar mesmo assim?":
        "'All' is not recommended. Continue anyway?",
    "Quantos resíduos? (1=Top 5, 2=Top 10)": "How many residues? (1=Top 5, 2=Top 10)",
    "  Ligante detectado em cada sistema:": "  Ligand detected in each system:",
    "não detectado": "not detected",
    "Manter o ligante no ranking? (mantido, ele tende a dominar o topo)":
        "Keep the ligand in the ranking? (if kept, it tends to dominate the top)",
    # cores
    "  Hex separados por vírgula (ex.: #0072B2,#D55E00). Enter usa a paleta padrão.":
        "  Comma-separated hex (e.g. #0072B2,#D55E00). Enter uses the default palette.",
    "Cores dos sistemas, na ordem escolhida": "System colors, in the chosen order",
    "Cores dos 5 componentes (vDW, elet., polar, apolar, total)":
        "Colors of the 5 components (vDW, elec., polar, nonpolar, total)",
    "Menos cores do que sistemas; o restante usa a paleta padrão.":
        "Fewer colors than systems; the rest use the default palette.",
    "Menos de 5 cores de componentes; completando com o padrão.":
        "Fewer than 5 component colors; filling in with the defaults.",
    # saída
    "  Organização dos arquivos gerados:": "  Organization of generated files:",
    "por tipo de gráfico (recomendado)": "by plot type (recommended)",
    "por formato": "by file format",
    "tudo na mesma pasta": "everything in one folder",
    "Organização": "Organization",
    "Formatos (png, pdf, svg; Enter = todos)": "Formats (png, pdf, svg; Enter = all)",
    "Exportar também os dados de energia em CSV (ns e kcal/mol)?":
        "Also export the energy data as CSV (ns and kcal/mol)?",
    "Abrir janela interativa para ajustar cada figura antes de salvar?":
        "Open an interactive window to adjust each figure before saving?",
    # geração
    "Gerando": "Generating",
    "Energia no tempo": "Energy vs. time",
    "Médias dos componentes": "Component averages",
    "Decomposição por resíduo": "Per-residue decomposition",
    "{label}: sem {fname}; ignorado nos gráficos de tempo.":
        "{label}: no {fname}; skipped in the time plots.",
    "{label}: sem {fname}; ignorado nas médias.":
        "{label}: no {fname}; skipped in the averages.",
    "{path}: nenhum frame da decomposição na faixa {a:g}–{b:g} ns":
        "{path}: no decomposition frames in the {a:g}–{b:g} ns range",
    "{label}: {fname} já vem agregado (média ± DP); resíduos calculados na trajetória inteira, não na faixa {a:g}–{b:g} ns.":
        "{label}: {fname} is already aggregated (mean ± SD); residues computed over the whole trajectory, not the {a:g}–{b:g} ns range.",
    "{label}: sem {fname}; ignorado nos resíduos.":
        "{label}: no {fname}; skipped in the residues.",
    "{label}: nenhum frame entre {a} e {b} ns (dados: {lo:.1f}–{hi:.1f} ns).":
        "{label}: no frames between {a} and {b} ns (data: {lo:.1f}–{hi:.1f} ns).",
    "Gráfico '{key}' não gerado: sem dados na faixa.":
        "Plot '{key}' not generated: no data in the range.",
    "Sem backend gráfico disponível; janela interativa ignorada.":
        "No graphical backend available; interactive window skipped.",
    "  ... ajuste '{stem}' na janela e feche-a para salvar.":
        "  ... adjust '{stem}' in the window and close it to save.",
    # início
    "gráficos de MM/PBSA e MM/GBSA (gmx_MMPBSA / Amber)": "MM/PBSA and MM/GBSA plots (gmx_MMPBSA / Amber)",
    "  Pasta raiz: {root}": "  Root folder: {root}",
    "Pasta não encontrada: {root}": "Folder not found: {root}",
    "Nenhum arquivo ({a}, {b}, {c}) encontrado em {root} ou nas subpastas.":
        "No file ({a}, {b}, {c}) found in {root} or its subfolders.",
    "  Fonte '{fam}' não instalada; usando '{used}'.":
        "  Font '{fam}' not installed; using '{used}'.",
    "  Configuração antiga atualizada (fonte Arial, textos por idioma).":
        "  Old configuration updated (Arial font, per-language texts).",
    # fim
    "Interrompido.": "Interrupted.",
    "INTERROMPIDO pelo usuário": "INTERRUPTED by the user",
    "Erro: {err}": "Error: {err}",
    "ERRO: {err}": "ERROR: {err}",
    "ERRO inesperado (veja o terminal)": "UNEXPECTED ERROR (see the terminal)",
    "Concluído.": "Done.",
    "{f} arquivo(s) de figura, {t} tabela(s)": "{f} figure file(s), {t} table(s)",
    "{w} aviso(s)": "{w} warning(s)",
    "Saída:": "Output:",
    "Relatório:": "Report:",
    # entrada
    "Argumento inválido para '{text}': {err}": "Invalid argument for '{text}': {err}",
    "Modo --yes: '{text}' não tem padrão; informe por argumento.":
        "--yes mode: '{text}' has no default; provide it as an argument.",
    "responda s ou n": "answer y or n",
    "'{part}' não é um número válido": "'{part}' is not a valid number",
    "{k} fora da faixa 1-{n}": "{k} is out of range 1-{n}",
    "nada selecionado": "nothing selected",
    "'{part}' não é um hex válido (ex.: #1f77b4)": "'{part}' is not a valid hex (e.g. #1f77b4)",
    "formato '{part}' não suportado (png, pdf, svg)":
        "format '{part}' not supported (png, pdf, svg)",
    "informe dois números: início e fim (ex.: 20 100)":
        "enter two numbers: start and end (e.g. 20 100)",
    "o início precisa ser menor que o fim": "the start must be smaller than the end",
    "escolha 1 (Top 5), 2 (Top 10) ou 3 (Todos)": "choose 1 (Top 5), 2 (Top 10) or 3 (All)",
    "escolha 1 (Top 5) ou 2 (Top 10)": "choose 1 (Top 5) or 2 (Top 10)",
    "escolha 1 (por tipo), 2 (por formato) ou 3 (tudo junto)":
        "choose 1 (by type), 2 (by format) or 3 (all together)",
    "nome vazio": "empty name",
    # erros de dados
    "{path} está malformado: {err}": "{path} is malformed: {err}",
    "{path}: não consegui ler os números ({err})": "{path}: could not read the numbers ({err})",
    "{path}: arquivo sem dados": "{path}: file has no data",
    "{path} tem {n} colunas; esperado: tempo + 5 componentes (vdW, elet., polar, apolar, total)":
        "{path} has {n} columns; expected: time + 5 components (vdW, elec., polar, nonpolar, total)",
    "{path}: faltam as colunas {missing}": "{path}: missing columns {missing}",
    "{path}: esperadas 3 colunas (nome, média, desvio)":
        "{path}: expected 3 columns (name, mean, deviation)",
    "{path}: valor não numérico na linha '{row}'": "{path}: non-numeric value in row '{row}'",
    "{path}: nenhum componente reconhecido": "{path}: no recognized component",
    # perguntas por prompt (ask)
    "Idioma / Language": "Idioma / Language",
}


EN_UI.update({
    "Energia no tempo — Fase gasosa (GGAS)": "Energy vs. time — Gas phase (GGAS)",
    "Energia no tempo — Solvatação (GSOLV)": "Energy vs. time — Solvation (GSOLV)",
    "  ⚠ {name}: ΔG total difere em {gap:.2f} kcal/mol (média) da soma vdW+elet.+polar+apolar "
    "(termos de ligação/1-4 não nulos no delta?).":
        "  ⚠ {name}: total ΔG differs by {gap:.2f} kcal/mol (mean) from vdW+elec.+polar+nonpolar "
        "(non-zero bonded/1-4 terms in the delta?).",
    "Modelo": "Model",
    "Qual método foi usado? (1=MM/GBSA, 2=MM/PBSA)": "Which method was used? (1=MM/GBSA, 2=MM/PBSA)",
    "escolha 1 (MM/GBSA) ou 2 (MM/PBSA)": "choose 1 (MM/GBSA) or 2 (MM/PBSA)",
    "não encontrado nos arquivos": "not found in the files",
    "Os arquivos não têm resultados {name}; encontrado: {found}.":
        "The files have no {name} results; found: {found}.",
    "Tempo em ns de cada frame da trajetória (coluna 'Frame #' do CSV)":
        "Time in ns of each trajectory frame ('Frame #' column of the CSV)",
    "informe um número, ex.: 0.01": "enter a number, e.g. 0.01",
    "deve ser maior que zero": "must be greater than zero",
    "{path}: não encontrei a tabela 'Delta Energy Terms' (rode o gmx_MMPBSA com -eo)":
        "{path}: 'Delta Energy Terms' table not found (run gmx_MMPBSA with -eo)",
    "{path}: não encontrei 'Total Energy Decomposition' na seção DELTAS (rode com -deo e &decomp)":
        "{path}: 'Total Energy Decomposition' not found in the DELTAS section (run with -deo and &decomp)",
})


if __name__ == "__main__":
    sys.exit(main())
