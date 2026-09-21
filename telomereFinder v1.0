#!/usr/bin/env python3
import re
import os
import csv
import threading
import queue
from collections import defaultdict
from datetime import datetime

import tkinter as tk
from tkinter import ttk, filedialog, messagebox


ORANGE_DARK = "#E65100"
ORANGE_MAIN = "#FB8C00"
ORANGE_MID = "#FFA726"
ORANGE_LIGHT = "#FFE0B2"
ORANGE_PALE = "#FFF3E0"
CREAM = "#FFFBF5"
BROWN_TEXT = "#3E2723"
BROWN_LIGHT = "#6D4C41"
WHITE = "#FFFFFF"
GREEN_OK = "#2E7D32"
RED_ERR = "#C62828"


TELOMERE_PRESETS = {
    'Plants': {
        'forward': ['CCCTAAA', 'CCCTAA'],
        'reverse': ['TTTAGGG', 'TTAGGG'],
        'description': 'Arabidopsis-type plant telomeres (CCCTAAA / TTTAGGG)',
    },
    'Vertebrates': {
        'forward': ['TTAGGG'],
        'reverse': ['CCCTAA'],
        'description': 'Vertebrate telomeres (TTAGGG / CCCTAA) - mammals, birds, fish, reptiles, amphibians',
    },
    'Insects (non-Drosophila)': {
        'forward': ['TTAGG'],
        'reverse': ['CCTAA'],
        'description': 'Insect telomeres (TTAGG / CCTAA) - excluding Drosophila, which uses retrotransposons',
    },
    'Nematodes': {
        'forward': ['TTAGGC'],
        'reverse': ['GCCTAA'],
        'description': 'Nematode telomeres (TTAGGC / GCCTAA) - includes Caenorhabditis elegans',
    },
    'Ciliates': {
        'forward': ['TTGGGG'],
        'reverse': ['CCCCAA'],
        'description': 'Ciliate telomeres (TTGGGG / CCCCAA) - includes Tetrahymena',
    },
    'Fungi (Candida)': {
        'forward': ['TGTGGG'],
        'reverse': ['CCCACA'],
        'description': 'Candida-type fungal telomeres (TGTGGG / CCCACA)',
    },
    'Fungi (Neurospora)': {
        'forward': ['TTAGGG'],
        'reverse': ['CCCTAA'],
        'description': 'Neurospora-type fungal telomeres (TTAGGG / CCCTAA)',
    },
    'Trypanosomes': {
        'forward': ['TTAGGG'],
        'reverse': ['CCCTAA'],
        'description': 'Trypanosome telomeres (TTAGGG / CCCTAA)',
    },
    'Custom': {
        'forward': [],
        'reverse': [],
        'description': 'User-defined telomeric motifs - enter comma-separated ACGT sequences',
    },
}


IUPAC = {
    'A': '[AG]',
    'G': '[AG]',
    'C': '[CT]',
    'T': '[CT]',
}


def generate_fuzzy_variants(motif):
    if not re.fullmatch(r'[ACGT]+', motif):
        return []
    if len(motif) < 5:
        return []
    positions = [1, 4]
    variants = []
    for pos in positions:
        if pos >= len(motif):
            continue
        base = motif[pos]
        if base in IUPAC:
            variants.append(motif[:pos] + IUPAC[base] + motif[pos + 1:])
    return variants


def parse_motif_list(text):
    if not text:
        return []
    parts = re.split(r'[,\s;]+', text.strip())
    return [p.upper() for p in parts if p]


def validate_motif(motif):
    return bool(re.fullmatch(r'[ACGTN]+', motif))


def suggest_preset(sequences, margin=2000, min_repeats=3):
    terminal_regions = []
    for seq in sequences.values():
        if len(seq) < 2 * margin:
            terminal_regions.append(seq)
        else:
            terminal_regions.append(seq[:margin])
            terminal_regions.append(seq[-margin:])

    scores = {}
    for preset_name, preset in TELOMERE_PRESETS.items():
        if preset_name == 'Custom':
            continue
        motifs = preset['forward'] + preset['reverse']
        if not motifs:
            continue
        total = 0
        for motif in motifs:
            try:
                regex = re.compile(f'({motif}){{{min_repeats},}}', re.IGNORECASE)
            except re.error:
                continue
            for region in terminal_regions:
                total += len(regex.findall(region))
        scores[preset_name] = total

    if not scores:
        return None, {}
    best_name = max(scores, key=scores.get)
    if scores[best_name] == 0:
        return None, scores
    return best_name, scores


class TelomereAnalyzer:
    def __init__(self, forward_motifs, reverse_motifs, min_repeats=3,
                 min_its_length=50, margin=5000, max_its=None, use_fuzzy=True):
        self.forward_motifs = [m.upper() for m in forward_motifs if m]
        self.reverse_motifs = [m.upper() for m in reverse_motifs if m]
        self.min_repeats = min_repeats
        self.min_its_length = min_its_length
        self.margin = margin
        self.max_its = max_its
        self.use_fuzzy = use_fuzzy
        self.patterns = self._build_patterns()

    def _build_patterns(self):
        patterns = []
        for motif in self.forward_motifs:
            patterns.append((motif, motif, 'F', False))
            if self.use_fuzzy:
                for v in generate_fuzzy_variants(motif):
                    patterns.append((v, motif + '_fuzzy', 'F', True))
        for motif in self.reverse_motifs:
            patterns.append((motif, motif, 'R', False))
            if self.use_fuzzy:
                for v in generate_fuzzy_variants(motif):
                    patterns.append((v, motif + '_fuzzy', 'R', True))
        return patterns

    def find_regions(self, sequence):
        regions = []
        for regex_str, name, orientation, is_fuzzy in self.patterns:
            try:
                regex = re.compile(f'({regex_str}){{{self.min_repeats},}}', re.IGNORECASE)
            except re.error:
                continue
            for match in regex.finditer(sequence):
                start, end = match.start(), match.end()
                length = end - start
                if length < self.min_its_length:
                    continue
                regions.append({
                    'start': start,
                    'end': end,
                    'length': length,
                    'seq': match.group(),
                    'orientation': orientation,
                    'pattern': name,
                    'is_fuzzy': is_fuzzy,
                })
        return self._deduplicate(regions)

    def _deduplicate(self, regions):
        regions = sorted(regions, key=lambda r: (r['start'], -r['length']))
        result = []
        for r in regions:
            merged = False
            for i, ex in enumerate(result):
                if r['start'] < ex['end'] and r['end'] > ex['start']:
                    if r['length'] > ex['length']:
                        result[i] = r
                    merged = True
                    break
            if not merged:
                result.append(r)
        return sorted(result, key=lambda r: r['start'])

    def classify(self, sequence):
        regions = self.find_regions(sequence)
        seq_len = len(sequence)
        start_ter = []
        end_ter = []
        internal = []
        for r in regions:
            if r['start'] <= self.margin:
                start_ter.append(r)
            elif (seq_len - r['end']) <= self.margin:
                end_ter.append(r)
            else:
                internal.append(r)
        return start_ter, end_ter, internal

    def reorganize(self, sequence, internal_regions):
        forward = sorted([r for r in internal_regions if r['orientation'] == 'F'],
                         key=lambda r: r['start'])
        reverse = sorted([r for r in internal_regions if r['orientation'] == 'R'],
                         key=lambda r: r['start'])

        if self.max_its:
            forward = sorted(sorted(forward, key=lambda r: -r['length'])[:self.max_its],
                             key=lambda r: r['start'])
            reverse = sorted(sorted(reverse, key=lambda r: -r['length'])[:self.max_its],
                             key=lambda r: r['start'])

        to_remove = sorted(forward + reverse, key=lambda r: r['start'], reverse=True)

        body = sequence
        for r in to_remove:
            body = body[:r['start']] + body[r['end']:]

        new_seq = ''.join(r['seq'] for r in forward) + body + ''.join(r['seq'] for r in reverse)

        log = []
        for r in forward:
            log.append({**r, 'destination': 'START'})
        for r in reverse:
            log.append({**r, 'destination': 'END'})

        return new_seq, log


def read_fasta(path):
    sequences = {}
    order = []
    current_id = None
    current_seq = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                order.append(current_id)
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    return sequences, order


def write_fasta(path, sequences, order, log_map):
    with open(path, 'w') as f:
        for sid in order:
            seq = sequences[sid]
            entries = log_map.get(sid, [])
            if entries:
                fwd = sum(1 for e in entries if e['destination'] == 'START')
                rev = sum(1 for e in entries if e['destination'] == 'END')
                f.write(f">{sid} [REORGANIZED F:{fwd} R:{rev}]\n")
            else:
                f.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def build_report(input_path, output_path, params, scaffold_stats, all_moves, elapsed):
    lines = []
    lines.append("=" * 80)
    lines.append("TELOMERE SCAFFOLD REORGANIZATION REPORT")
    lines.append("=" * 80)
    lines.append(f"Generated  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Input      : {input_path}")
    lines.append(f"Output     : {output_path}")
    lines.append(f"Elapsed    : {elapsed:.2f} s")
    lines.append("")
    lines.append("ORGANISM AND MOTIFS")
    lines.append("-" * 80)
    lines.append(f"Organism preset      : {params.get('preset', 'Custom')}")
    lines.append(f"Forward motifs       : {', '.join(params.get('forward_motifs', []))}")
    lines.append(f"Reverse motifs       : {', '.join(params.get('reverse_motifs', []))}")
    lines.append(f"Fuzzy variants       : {'enabled' if params.get('use_fuzzy') else 'disabled'}")
    lines.append("")
    lines.append("DETECTION PARAMETERS")
    lines.append("-" * 80)
    lines.append(f"Min repeats per array        : {params['min_repeats']}")
    lines.append(f"Min ITS length (bp)          : {params['min_its_length']}")
    lines.append(f"Terminal margin (bp)         : {params['margin']}")
    lines.append(f"Max ITS per scaffold         : {params['max_its'] if params['max_its'] else 'unlimited'}")
    lines.append("")

    total_scaffolds = len(scaffold_stats)
    with_its = sum(1 for s in scaffold_stats.values() if s['moved'] > 0)
    total_moved = sum(s['moved'] for s in scaffold_stats.values())
    total_fwd = sum(s['forward'] for s in scaffold_stats.values())
    total_rev = sum(s['reverse'] for s in scaffold_stats.values())
    total_bp = sum(s['bp_moved'] for s in scaffold_stats.values())
    total_start_ter = sum(s['start_ter'] for s in scaffold_stats.values())
    total_end_ter = sum(s['end_ter'] for s in scaffold_stats.values())

    lines.append("GLOBAL SUMMARY")
    lines.append("-" * 80)
    lines.append(f"Total scaffolds analyzed          : {total_scaffolds}")
    lines.append(f"Scaffolds with terminal telomeres : {total_start_ter} start / {total_end_ter} end")
    lines.append(f"Scaffolds with internal telomeres : {with_its}")
    lines.append(f"Total ITSs relocated              : {total_moved}")
    lines.append(f"  Forward moved to START          : {total_fwd}")
    lines.append(f"  Reverse moved to END            : {total_rev}")
    lines.append(f"Total bases relocated             : {total_bp:,} bp")
    lines.append("")

    lines.append("PER-SCAFFOLD SUMMARY")
    lines.append("-" * 80)
    header = f"{'Scaffold':<24} {'Length':>12} {'Start_ter':>10} {'End_ter':>8} {'ITS_F':>6} {'ITS_R':>6} {'Bases':>10}"
    lines.append(header)
    lines.append("-" * 80)
    for sid in scaffold_stats:
        s = scaffold_stats[sid]
        lines.append(
            f"{sid:<24} {s['length']:>12,} {s['start_ter']:>10} {s['end_ter']:>8} "
            f"{s['forward']:>6} {s['reverse']:>6} {s['bp_moved']:>10,}"
        )
    lines.append("")

    if all_moves:
        lines.append("DETAILED ITS MOVEMENTS")
        lines.append("-" * 80)
        lines.append(f"{'Scaffold':<24} {'Start':>10} {'End':>10} {'Length':>8} {'Ori':>5} {'Dest':>7} {'Pattern':<18}")
        lines.append("-" * 80)
        for m in all_moves:
            lines.append(
                f"{m['scaffold']:<24} {m['start']:>10,} {m['end']:>10,} "
                f"{m['length']:>8} {m['orientation']:>5} {m['destination']:>7} {m['pattern']:<18}"
            )
        lines.append("")

    lines.append("END OF REPORT")
    lines.append("=" * 80)

    return "\n".join(lines)


class TelomereReorganizerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Telomere Scaffold Reorganizer")
        self.root.geometry("1180x940")
        self.root.minsize(1040, 800)
        self.root.configure(bg=ORANGE_PALE)

        self.msg_queue = queue.Queue()
        self.worker = None
        self.suggest_worker = None
        self.running = False

        self.var_input = tk.StringVar()
        self.var_output = tk.StringVar()
        self.var_report = tk.StringVar()
        self.var_preset = tk.StringVar(value='Plants')
        self.var_use_fuzzy = tk.BooleanVar(value=True)
        self.var_forward_motifs = tk.StringVar()
        self.var_reverse_motifs = tk.StringVar()
        self.var_min_repeats = tk.StringVar(value="3")
        self.var_min_its = tk.StringVar(value="50")
        self.var_margin = tk.StringVar(value="5000")
        self.var_max_its = tk.StringVar(value="")
        self.var_status = tk.StringVar(value="Ready")
        self.var_preset_desc = tk.StringVar()

        self._setup_styles()
        self._build_ui()
        self._on_preset_change()
        self._poll_queue()

    def _setup_styles(self):
        style = ttk.Style()
        try:
            style.theme_use('clam')
        except tk.TclError:
            pass

        style.configure('TFrame', background=ORANGE_PALE)
        style.configure('Card.TFrame', background=CREAM)
        style.configure('Header.TFrame', background=ORANGE_DARK)
        style.configure('Status.TFrame', background=ORANGE_LIGHT)

        style.configure('TLabel', background=ORANGE_PALE, foreground=BROWN_TEXT,
                        font=('Segoe UI', 10))
        style.configure('Card.TLabel', background=CREAM, foreground=BROWN_TEXT,
                        font=('Segoe UI', 10))
        style.configure('HeaderTitle.TLabel', background=ORANGE_DARK, foreground=WHITE,
                        font=('Segoe UI', 20, 'bold'))
        style.configure('HeaderSub.TLabel', background=ORANGE_DARK, foreground=ORANGE_LIGHT,
                        font=('Segoe UI', 10))
        style.configure('SectionTitle.TLabel', background=ORANGE_PALE,
                        foreground=ORANGE_DARK, font=('Segoe UI', 11, 'bold'))
        style.configure('FieldLabel.TLabel', background=CREAM, foreground=BROWN_LIGHT,
                        font=('Segoe UI', 9))
        style.configure('Desc.TLabel', background=CREAM, foreground=BROWN_LIGHT,
                        font=('Segoe UI', 9, 'italic'))
        style.configure('Status.TLabel', background=ORANGE_LIGHT, foreground=BROWN_TEXT,
                        font=('Segoe UI', 9))

        style.configure('TEntry', fieldbackground=WHITE, foreground=BROWN_TEXT,
                        bordercolor=ORANGE_LIGHT, lightcolor=ORANGE_MID,
                        darkcolor=ORANGE_MID, insertcolor=ORANGE_DARK, padding=6)
        style.map('TEntry',
                  bordercolor=[('focus', ORANGE_MAIN)],
                  lightcolor=[('focus', ORANGE_MAIN)],
                  darkcolor=[('focus', ORANGE_MAIN)])

        style.configure('TCombobox', fieldbackground=WHITE, background=WHITE,
                        foreground=BROWN_TEXT, arrowcolor=ORANGE_DARK,
                        bordercolor=ORANGE_LIGHT, lightcolor=ORANGE_MID,
                        darkcolor=ORANGE_MID, padding=4)
        style.map('TCombobox',
                  fieldbackground=[('readonly', WHITE)],
                  bordercolor=[('focus', ORANGE_MAIN)])

        style.configure('TCheckbutton', background=CREAM, foreground=BROWN_TEXT,
                        font=('Segoe UI', 9))
        style.map('TCheckbutton', background=[('active', CREAM)])

        style.configure('Orange.TButton', background=ORANGE_MAIN, foreground=WHITE,
                        font=('Segoe UI', 10, 'bold'), padding=(18, 10),
                        borderwidth=0, relief='flat')
        style.map('Orange.TButton',
                  background=[('active', ORANGE_DARK), ('pressed', ORANGE_DARK),
                              ('disabled', ORANGE_LIGHT)],
                  foreground=[('disabled', '#AAAAAA')])

        style.configure('Outline.TButton', background=CREAM, foreground=ORANGE_DARK,
                        font=('Segoe UI', 10, 'bold'), padding=(16, 10),
                        borderwidth=1, relief='solid')
        style.map('Outline.TButton',
                  background=[('active', ORANGE_LIGHT)],
                  bordercolor=[('active', ORANGE_MAIN)])

        style.configure('Small.TButton', background=ORANGE_MID, foreground=WHITE,
                        font=('Segoe UI', 9, 'bold'), padding=(10, 6), borderwidth=0)
        style.map('Small.TButton',
                  background=[('active', ORANGE_DARK)])

        style.configure('Orange.Horizontal.TProgressbar',
                        background=ORANGE_MAIN, troughcolor=ORANGE_LIGHT,
                        bordercolor=ORANGE_LIGHT, lightcolor=ORANGE_MAIN,
                        darkcolor=ORANGE_MAIN, thickness=8)

    def _build_ui(self):
        header = ttk.Frame(self.root, style='Header.TFrame')
        header.pack(fill='x')

        header_inner = ttk.Frame(header, style='Header.TFrame')
        header_inner.pack(fill='x', padx=24, pady=16)

        ttk.Label(header_inner, text="\u25C9  telomereFinder",
                  style='HeaderTitle.TLabel').pack(anchor='w')
        ttk.Label(header_inner,
                  text="Detect and organize telomeric sequences in de novo assembled genomes",
                  style='HeaderSub.TLabel').pack(anchor='w', pady=(4, 0))

        body = ttk.Frame(self.root, style='TFrame')
        body.pack(fill='both', expand=True, padx=20, pady=16)

        files_card = ttk.Frame(body, style='Card.TFrame')
        files_card.pack(fill='x', pady=(0, 12))

        files_inner = ttk.Frame(files_card, style='Card.TFrame')
        files_inner.pack(fill='x', padx=18, pady=16)

        ttk.Label(files_inner, text="FILES", style='SectionTitle.TLabel',
                  background=CREAM, foreground=ORANGE_DARK).grid(
            row=0, column=0, columnspan=3, sticky='w', pady=(0, 10))

        self._file_row(files_inner, 1, "Input FASTA", self.var_input,
                       self._browse_input)
        self._file_row(files_inner, 2, "Output FASTA", self.var_output,
                       self._browse_output)
        self._file_row(files_inner, 3, "Report file", self.var_report,
                       self._browse_report)

        files_inner.columnconfigure(1, weight=1)

        org_card = ttk.Frame(body, style='Card.TFrame')
        org_card.pack(fill='x', pady=(0, 12))

        org_inner = ttk.Frame(org_card, style='Card.TFrame')
        org_inner.pack(fill='x', padx=18, pady=16)

        ttk.Label(org_inner, text="ORGANISM PRESET", style='SectionTitle.TLabel',
                  background=CREAM, foreground=ORANGE_DARK).grid(
            row=0, column=0, columnspan=4, sticky='w', pady=(0, 10))

        ttk.Label(org_inner, text="Organism group", style='FieldLabel.TLabel').grid(
            row=1, column=0, sticky='w', pady=5, padx=(0, 12))

        self.preset_combo = ttk.Combobox(org_inner, textvariable=self.var_preset,
                                          values=list(TELOMERE_PRESETS.keys()),
                                          state='readonly', width=32)
        self.preset_combo.grid(row=1, column=1, sticky='w', pady=5, padx=(0, 8))
        self.preset_combo.bind('<<ComboboxSelected>>', self._on_preset_change)

        self.btn_suggest = ttk.Button(org_inner, text="Suggest Motifs",
                                       style='Small.TButton',
                                       command=self._suggest_motifs)
        self.btn_suggest.grid(row=1, column=2, sticky='w', pady=5, padx=(0, 8))

        self.var_preset_desc.set(TELOMERE_PRESETS['Plants']['description'])
        ttk.Label(org_inner, textvariable=self.var_preset_desc, style='Desc.TLabel',
                  wraplength=680).grid(row=2, column=0, columnspan=4,
                                        sticky='w', pady=(0, 10))

        ttk.Label(org_inner, text="Forward motifs", style='FieldLabel.TLabel').grid(
            row=3, column=0, sticky='w', pady=5, padx=(0, 12))
        self.entry_fwd = ttk.Entry(org_inner, textvariable=self.var_forward_motifs)
        self.entry_fwd.grid(row=3, column=1, columnspan=3, sticky='ew',
                            pady=5, padx=(0, 8))

        ttk.Label(org_inner, text="Reverse motifs", style='FieldLabel.TLabel').grid(
            row=4, column=0, sticky='w', pady=5, padx=(0, 12))
        self.entry_rev = ttk.Entry(org_inner, textvariable=self.var_reverse_motifs)
        self.entry_rev.grid(row=4, column=1, columnspan=3, sticky='ew',
                            pady=5, padx=(0, 8))

        ttk.Checkbutton(org_inner, text="Use fuzzy variants (auto-generated from canonical motifs)",
                        variable=self.var_use_fuzzy).grid(
            row=5, column=0, columnspan=4, sticky='w', pady=(8, 0))

        org_inner.columnconfigure(1, weight=1)
        org_inner.columnconfigure(3, weight=1)

        params_card = ttk.Frame(body, style='Card.TFrame')
        params_card.pack(fill='x', pady=(0, 12))

        params_inner = ttk.Frame(params_card, style='Card.TFrame')
        params_inner.pack(fill='x', padx=18, pady=16)

        ttk.Label(params_inner, text="DETECTION PARAMETERS", style='SectionTitle.TLabel',
                  background=CREAM, foreground=ORANGE_DARK).grid(
            row=0, column=0, columnspan=4, sticky='w', pady=(0, 10))

        self._param_field(params_inner, 1, 0, "Min repeats", self.var_min_repeats)
        self._param_field(params_inner, 1, 2, "Min ITS length (bp)", self.var_min_its)
        self._param_field(params_inner, 2, 0, "Terminal margin (bp)", self.var_margin)
        self._param_field(params_inner, 2, 2, "Max ITS per scaffold", self.var_max_its,
                          hint="empty = unlimited")

        params_inner.columnconfigure(1, weight=1)
        params_inner.columnconfigure(3, weight=1)

        actions = ttk.Frame(body, style='TFrame')
        actions.pack(fill='x', pady=(0, 12))

        self.btn_run = ttk.Button(actions, text="Reorganize", style='Orange.TButton',
                                  command=lambda: self._start(dry_run=False))
        self.btn_run.pack(side='left', padx=(0, 8))

        self.btn_dry = ttk.Button(actions, text="Dry Run", style='Outline.TButton',
                                  command=lambda: self._start(dry_run=True))
        self.btn_dry.pack(side='left', padx=(0, 8))

        self.btn_clear = ttk.Button(actions, text="Clear Log", style='Outline.TButton',
                                    command=self._clear_log)
        self.btn_clear.pack(side='left')

        self.progress = ttk.Progressbar(body, style='Orange.Horizontal.TProgressbar',
                                        mode='determinate', maximum=100)
        self.progress.pack(fill='x', pady=(0, 12))

        log_card = ttk.Frame(body, style='Card.TFrame')
        log_card.pack(fill='both', expand=True)

        log_inner = ttk.Frame(log_card, style='Card.TFrame')
        log_inner.pack(fill='both', expand=True, padx=18, pady=16)

        ttk.Label(log_inner, text="ANALYSIS LOG", style='SectionTitle.TLabel',
                  background=CREAM, foreground=ORANGE_DARK).pack(anchor='w', pady=(0, 8))

        log_wrap = tk.Frame(log_inner, bg=ORANGE_LIGHT, bd=0,
                            highlightthickness=1, highlightbackground=ORANGE_LIGHT)
        log_wrap.pack(fill='both', expand=True)

        self.log_text = tk.Text(log_wrap, wrap='word', bg=WHITE, fg=BROWN_TEXT,
                                font=('Consolas', 9), relief='flat', bd=0,
                                padx=10, pady=8, insertbackground=ORANGE_DARK)
        self.log_text.pack(side='left', fill='both', expand=True)

        scroll = ttk.Scrollbar(log_wrap, orient='vertical',
                               command=self.log_text.yview)
        scroll.pack(side='right', fill='y')
        self.log_text.configure(yscrollcommand=scroll.set)

        self.log_text.tag_configure('info', foreground=BROWN_TEXT)
        self.log_text.tag_configure('ok', foreground=GREEN_OK)
        self.log_text.tag_configure('err', foreground=RED_ERR)
        self.log_text.tag_configure('head', foreground=ORANGE_DARK,
                                    font=('Consolas', 9, 'bold'))
        self.log_text.tag_configure('dim', foreground=BROWN_LIGHT)

        status_bar = ttk.Frame(self.root, style='Status.TFrame')
        status_bar.pack(fill='x', side='bottom')
        status_inner = ttk.Frame(status_bar, style='Status.TFrame')
        status_inner.pack(fill='x', padx=20, pady=6)
        ttk.Label(status_inner, textvariable=self.var_status,
                  style='Status.TLabel').pack(side='left')

    def _file_row(self, parent, row, label, var, browse_cmd):
        ttk.Label(parent, text=label, style='FieldLabel.TLabel').grid(
            row=row, column=0, sticky='w', pady=5, padx=(0, 12))
        entry = ttk.Entry(parent, textvariable=var)
        entry.grid(row=row, column=1, sticky='ew', pady=5, padx=(0, 8))
        ttk.Button(parent, text="Browse", style='Small.TButton',
                   command=browse_cmd).grid(row=row, column=2, pady=5)

    def _param_field(self, parent, row, col, label, var, hint=None):
        ttk.Label(parent, text=label, style='FieldLabel.TLabel').grid(
            row=row, column=col, sticky='w', pady=5, padx=(0, 12))
        frame = ttk.Frame(parent, style='Card.TFrame')
        frame.grid(row=row, column=col + 1, sticky='ew', pady=5, padx=(0, 20))
        entry = ttk.Entry(frame, textvariable=var, width=18)
        entry.pack(side='left', fill='x', expand=True)
        if hint:
            ttk.Label(frame, text=f"  {hint}", style='FieldLabel.TLabel').pack(side='left')

    def _on_preset_change(self, event=None):
        name = self.var_preset.get()
        preset = TELOMERE_PRESETS.get(name)
        if not preset:
            return
        self.var_preset_desc.set(preset['description'])
        if name == 'Custom':
            self.entry_fwd.configure(state='normal')
            self.entry_rev.configure(state='normal')
        else:
            self.var_forward_motifs.set(', '.join(preset['forward']))
            self.var_reverse_motifs.set(', '.join(preset['reverse']))
            self.entry_fwd.configure(state='readonly')
            self.entry_rev.configure(state='readonly')

    def _browse_input(self):
        path = filedialog.askopenfilename(
            title="Select input FASTA",
            filetypes=[("FASTA files", "*.fasta *.fa *.fas *.fna *.ffn"),
                       ("All files", "*.*")])
        if not path:
            return
        self.var_input.set(path)
        base, ext = os.path.splitext(path)
        self.var_output.set(f"{base}_reorganized{ext if ext else '.fasta'}")
        self.var_report.set(f"{base}_telomere_report.txt")

    def _browse_output(self):
        path = filedialog.asksaveasfilename(
            title="Save reorganized FASTA",
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta *.fa *.fas"), ("All files", "*.*")])
        if path:
            self.var_output.set(path)

    def _browse_report(self):
        path = filedialog.asksaveasfilename(
            title="Save report file",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if path:
            self.var_report.set(path)

    def _log(self, text, tag='info'):
        self.log_text.insert('end', text + "\n", tag)
        self.log_text.see('end')

    def _clear_log(self):
        self.log_text.delete('1.0', 'end')

    def _set_running(self, running):
        self.running = running
        state = 'disabled' if running else 'normal'
        self.btn_run.configure(state=state)
        self.btn_dry.configure(state=state)
        self.progress['value'] = 0

    def _queue_log(self, text, tag='info'):
        self.msg_queue.put(('log', text, tag))

    def _queue_progress(self, value):
        self.msg_queue.put(('progress', value))

    def _queue_status(self, text):
        self.msg_queue.put(('status', text))

    def _queue_done(self, success, message):
        self.msg_queue.put(('done', success, message))

    def _queue_suggest(self, best_name, scores):
        self.msg_queue.put(('suggest_done', best_name, scores))

    def _poll_queue(self):
        try:
            while True:
                msg = self.msg_queue.get_nowait()
                if msg[0] == 'log':
                    self._log(msg[1], msg[2])
                elif msg[0] == 'progress':
                    self.progress['value'] = msg[1]
                elif msg[0] == 'status':
                    self.var_status.set(msg[1])
                elif msg[0] == 'done':
                    self._set_running(False)
                    if msg[1]:
                        self.var_status.set("Completed")
                    else:
                        self.var_status.set("Failed")
                        messagebox.showerror("Error", msg[2])
                elif msg[0] == 'suggest_done':
                    self.btn_suggest.configure(state='normal')
                    self.var_status.set("Ready")
                    if msg[1]:
                        self.var_preset.set(msg[1])
                        self._on_preset_change()
                        self._log(f"Suggested preset: {msg[1]}", 'ok')
                        top = sorted(msg[2].items(), key=lambda x: -x[1])[:3]
                        for name, score in top:
                            self._log(f"  {name:<30} score={score}", 'dim')
                    else:
                        self._log("No preset could be confidently suggested from terminal regions.", 'err')
        except queue.Empty:
            pass
        self.root.after(80, self._poll_queue)

    def _suggest_motifs(self):
        input_path = self.var_input.get().strip()
        if not input_path or not os.path.isfile(input_path):
            messagebox.showwarning("Missing input", "Please select a valid input FASTA file first.")
            return
        if self.running:
            return
        try:
            min_repeats = int(self.var_min_repeats.get())
        except ValueError:
            min_repeats = 3
        self._log("")
        self._log("=" * 72, 'head')
        self._log("Scanning scaffold ends to suggest organism preset...", 'head')
        self.btn_suggest.configure(state='disabled')
        self.var_status.set("Scanning terminal regions...")
        self.suggest_worker = threading.Thread(
            target=self._run_suggest,
            args=(input_path, min_repeats),
            daemon=True
        )
        self.suggest_worker.start()

    def _run_suggest(self, input_path, min_repeats):
        try:
            self._queue_log("Reading FASTA for motif suggestion...", 'head')
            sequences, _ = read_fasta(input_path)
            self._queue_log(f"  Loaded {len(sequences)} sequences", 'ok')
            self._queue_log(f"  Scanning terminal regions (margin=2000 bp, min_repeats={min_repeats})...")
            best, scores = suggest_preset(sequences, margin=2000, min_repeats=min_repeats)
            if best:
                self._queue_log(f"Best match: {best}", 'ok')
            else:
                self._queue_log("No tandem arrays found at scaffold ends.", 'err')
            self._queue_suggest(best, scores)
        except Exception as exc:
            self._queue_log(f"ERROR during suggestion: {exc}", 'err')
            self._queue_suggest(None, {})

    def _start(self, dry_run):
        if self.running:
            return

        input_path = self.var_input.get().strip()
        if not input_path or not os.path.isfile(input_path):
            messagebox.showwarning("Missing input", "Please select a valid input FASTA file.")
            return

        output_path = self.var_output.get().strip()
        report_path = self.var_report.get().strip()

        if not dry_run:
            if not output_path:
                messagebox.showwarning("Missing output", "Please specify an output FASTA path.")
                return
            if not report_path:
                messagebox.showwarning("Missing report", "Please specify a report file path.")
                return

        preset_name = self.var_preset.get()
        forward_motifs = parse_motif_list(self.var_forward_motifs.get())
        reverse_motifs = parse_motif_list(self.var_reverse_motifs.get())

        if not forward_motifs and not reverse_motifs:
            messagebox.showerror("Missing motifs",
                                 "Please provide at least one forward or reverse motif.")
            return

        for m in forward_motifs + reverse_motifs:
            if not validate_motif(m):
                messagebox.showerror("Invalid motif",
                                     f"Motif '{m}' contains invalid characters. Use only A, C, G, T, N.")
                return

        try:
            min_repeats = int(self.var_min_repeats.get())
            min_its = int(self.var_min_its.get())
            margin = int(self.var_margin.get())
        except ValueError:
            messagebox.showerror("Invalid parameters",
                                 "Min repeats, Min ITS and Margin must be integers.")
            return

        max_its_raw = self.var_max_its.get().strip()
        max_its = None
        if max_its_raw:
            try:
                max_its = int(max_its_raw)
                if max_its <= 0:
                    max_its = None
            except ValueError:
                messagebox.showerror("Invalid parameter", "Max ITS must be an integer.")
                return

        if min_repeats < 2:
            messagebox.showerror("Invalid parameter", "Min repeats must be >= 2.")
            return
        if min_its < 1:
            messagebox.showerror("Invalid parameter", "Min ITS length must be >= 1.")
            return
        if margin < 0:
            messagebox.showerror("Invalid parameter", "Margin must be >= 0.")
            return

        params = {
            'preset': preset_name,
            'forward_motifs': forward_motifs,
            'reverse_motifs': reverse_motifs,
            'use_fuzzy': self.var_use_fuzzy.get(),
            'min_repeats': min_repeats,
            'min_its_length': min_its,
            'margin': margin,
            'max_its': max_its,
        }

        self._set_running(True)
        self.var_status.set("Running...")
        self._log("=" * 72, 'head')
        self._log("TELOMERE SCAFFOLD REORGANIZATION", 'head')
        self._log("=" * 72, 'head')
        self._log(f"Input      : {input_path}")
        self._log(f"Mode       : {'DRY RUN' if dry_run else 'REORGANIZE'}")
        self._log(f"Preset     : {preset_name}")
        self._log(f"Fwd motifs : {', '.join(forward_motifs)}")
        self._log(f"Rev motifs : {', '.join(reverse_motifs)}")
        self._log(f"Fuzzy      : {'enabled' if params['use_fuzzy'] else 'disabled'}")
        self._log(f"Min repeats={min_repeats}  Min ITS={min_its}  "
                  f"Margin={margin}  Max ITS={max_its if max_its else 'unlimited'}")
        self._log("")

        self.worker = threading.Thread(
            target=self._run_pipeline,
            args=(input_path, output_path, report_path, params, dry_run),
            daemon=True
        )
        self.worker.start()

    def _run_pipeline(self, input_path, output_path, report_path, params, dry_run):
        start_time = datetime.now()
        try:
            self._queue_status("Reading FASTA...")
            self._queue_log("Reading input FASTA...", 'head')
            sequences, order = read_fasta(input_path)
            self._queue_log(f"  Loaded {len(sequences)} sequences", 'ok')

            analyzer = TelomereAnalyzer(
                forward_motifs=params['forward_motifs'],
                reverse_motifs=params['reverse_motifs'],
                min_repeats=params['min_repeats'],
                min_its_length=params['min_its_length'],
                margin=params['margin'],
                max_its=params['max_its'],
                use_fuzzy=params['use_fuzzy'],
            )

            self._queue_log(f"  Patterns loaded: {len(analyzer.patterns)}", 'dim')

            scaffold_stats = {}
            log_map = defaultdict(list)
            all_moves = []
            reorganized = {}
            total = len(order)
            moved_scaffolds = 0

            self._queue_log("")
            self._queue_log("Analyzing scaffolds and locating telomeric repeats...", 'head')

            for idx, sid in enumerate(order):
                seq = sequences[sid]
                start_ter, end_ter, internal = analyzer.classify(seq)

                stat = {
                    'length': len(seq),
                    'start_ter': len(start_ter),
                    'end_ter': len(end_ter),
                    'forward': 0,
                    'reverse': 0,
                    'moved': 0,
                    'bp_moved': 0,
                }

                if internal and not dry_run:
                    new_seq, moves = analyzer.reorganize(seq, internal)
                    reorganized[sid] = new_seq
                    for m in moves:
                        entry = dict(m)
                        entry['scaffold'] = sid
                        log_map[sid].append(entry)
                        all_moves.append(entry)
                        if m['destination'] == 'START':
                            stat['forward'] += 1
                        else:
                            stat['reverse'] += 1
                        stat['bp_moved'] += m['length']
                    stat['moved'] = len(moves)
                    if moves:
                        moved_scaffolds += 1
                    self._queue_log(
                        f"  {sid}: {len(start_ter)} start-ter, {len(end_ter)} end-ter, "
                        f"{len(internal)} ITS -> moved F:{stat['forward']} R:{stat['reverse']} "
                        f"({stat['bp_moved']:,} bp)", 'ok')
                elif internal and dry_run:
                    stat['moved'] = len(internal)
                    stat['bp_moved'] = sum(r['length'] for r in internal)
                    for m in internal:
                        entry = dict(m)
                        entry['scaffold'] = sid
                        entry['destination'] = 'START' if m['orientation'] == 'F' else 'END'
                        all_moves.append(entry)
                    moved_scaffolds += 1
                    self._queue_log(
                        f"  {sid}: {len(start_ter)} start-ter, {len(end_ter)} end-ter, "
                        f"{len(internal)} internal ITS (dry-run)", 'ok')
                else:
                    reorganized[sid] = seq

                scaffold_stats[sid] = stat
                self._queue_progress((idx + 1) / total * 100)
                self._queue_status(f"Analyzing... {idx + 1}/{total}")

            self._queue_log("")
            self._queue_log(
                f"Detected {len(all_moves)} internal telomere regions in "
                f"{moved_scaffolds} scaffolds", 'head')

            if dry_run:
                self._queue_log("Dry-run mode: no files were written", 'dim')
            else:
                self._queue_status("Writing reorganized FASTA...")
                self._queue_log("")
                self._queue_log(f"Writing reorganized FASTA: {output_path}", 'head')
                write_fasta(output_path, reorganized, order, log_map)
                self._queue_log("  Done", 'ok')

                self._queue_status("Writing report...")
                self._queue_log(f"Writing report: {report_path}", 'head')
                elapsed = (datetime.now() - start_time).total_seconds()
                report = build_report(input_path, output_path, params,
                                      scaffold_stats, all_moves, elapsed)
                with open(report_path, 'w') as f:
                    f.write(report)
                self._queue_log("  Done", 'ok')

            elapsed = (datetime.now() - start_time).total_seconds()
            self._queue_log("")
            self._queue_log("=" * 72, 'head')
            self._queue_log(f"Finished in {elapsed:.2f} s", 'head')
            self._queue_log(f"  Scaffolds total     : {len(order)}")
            self._queue_log(f"  Scaffolds affected  : {moved_scaffolds}")
            self._queue_log(f"  ITS regions found   : {len(all_moves)}")
            self._queue_log("=" * 72, 'head')

            self._queue_done(True, "")

        except Exception as exc:
            self._queue_log(f"ERROR: {exc}", 'err')
            self._queue_done(False, str(exc))


def main():
    root = tk.Tk()
    app = TelomereReorganizerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
