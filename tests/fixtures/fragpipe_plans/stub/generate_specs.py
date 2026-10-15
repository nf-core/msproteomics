"""Generate tests/fixtures/fragpipe_plans/stub/<name>.plan_inputs.json for every stub test reaching FRAGPIPE_TOOL_CONFIGS.

Run: python3 -I tests/fixtures/fragpipe_plans/stub/generate_specs.py <repo>; then record them with
tests/lib/record_fragpipe_plans.nf.test (README.md). A spec that no longer matches its test stops the stub test in
checkFragpipePlanManifest (file, experiment, bioreplicate, data type and FragPipe group are compared).

Mirrors parseSamplesheet / uniqueSpectraFiles / fragpipeGroup / guessFragpipeDataType / tmtAnnotationLines
(subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils.nf).
"""
import csv
import io
import json
import os
import re
import sys

REPO = sys.argv[1]
FIX = f'{REPO}/tests/fixtures'
OUT = f'{FIX}/fragpipe_plans/stub'

TMT16 = '126 127N 127C 128N 128C 129N 129C 130N 130C 131N 131C 132N 132C 133N 133C 134N'.split()
TMT_LABELS = {
    'TMT-6': '126, 127N, 128C, 129N, 130C, 131N'.split(', '),
    'TMT-10': '126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N'.split(', '),
}


def read(path):
    with open(path, encoding='utf-8', newline='') as fh:
        return fh.read()


def props(text):
    """java.util.Properties.load, enough for FragPipe .workflow files (no continuation lines)."""
    out = {}
    for line in text.splitlines():
        s = line.lstrip(' \t\f')
        if not s or s[0] in '#!':
            continue
        m = re.match(r'([^=:\s]+)\s*[=:]?\s*(.*)$', s)
        if m:
            out[m.group(1)] = m.group(2)
        else:
            out[s] = ''
    return out


def speclib_or_diann(p):
    return p.get('speclibgen.run-speclibgen') == 'true' or p.get('diann.run-dia-nn') == 'true'


def fp_experiment(e):
    return '' if e is None else re.sub(r'[^A-Za-z0-9_]', '_', e.strip())


def fp_group(e, b, single):
    if single:
        return ''
    exp = fp_experiment(e)
    if b is not None:
        return f'{exp}_{b}' if exp else f'exp_{b}'
    return exp


def guess_type(path):
    name = os.path.normpath(os.path.abspath(path))
    if 'dda' in name.lower() or any(t in name for t in ['_Q1.', '_Q2.', '_Q3.', '_diatracer.']):
        return 'DDA'
    if 'DIA' in name:
        return 'DIA'
    if 'WWA' in name:
        return 'DDA+'
    return 'DDA'


def parse_samplesheet(text, single):
    rows = []
    for cells in csv.DictReader(io.StringIO(text)):
        c = {k: (v if v != '' else None) for k, v in cells.items()}
        meta = {
            'sample': c['sample'],
            'condition': c.get('condition') if c.get('condition') is not None else c['sample'],
            'label': c.get('label'),
            'experiment': c.get('experiment'),
            'bioreplicate': int(c['bioreplicate']) if c.get('bioreplicate') is not None else None,
        }
        meta['group'] = fp_group(meta['experiment'], meta['bioreplicate'], single)
        meta['data_type'] = c.get('data_type') or guess_type(c['spectra'])
        rows.append((meta, c['spectra']))
    return rows


def unique_files(rows):
    seen = {}
    for meta, path in rows:
        seen.setdefault(path, (meta, path))
    return list(seen.values())


def annotation(rows, files, p):
    if p.get('tmtintegrator.run-tmtintegrator') != 'true':
        return []
    reagents = TMT_LABELS[p['tmtintegrator.channel_num'].strip()]
    out = []
    for g in sorted({m['group'] for m, _ in files}):
        labelled = [m for m, _ in rows if m['group'] == g and m['label']]
        if not labelled:
            out += [[g, ch, f'{g}_{ch}'] for ch in reagents]
            continue
        by = {}
        for m in labelled:
            by.setdefault(re.sub(r'^TMT', '', m['label']), m)
        assert sorted(by) == sorted(reagents), (g, by.keys())
        out += [[g, ch, by[ch]['sample']] for ch in reagents]
    return out


def rel(path):
    return os.path.relpath(path, REPO) if path.startswith(REPO) else path


def pipeline_spec(samplesheet_text, workflow_text, database_name, database_source, threads='--threads 2'):
    p = props(workflow_text)
    rows = parse_samplesheet(samplesheet_text, speclib_or_diann(p))
    files = unique_files(rows)
    return {
        'workflow_text': workflow_text,
        'manifest_rows': [[os.path.basename(path), m['experiment'], m['bioreplicate'], m['data_type'], m['group']] for m, path in files],
        'annotation_rows': annotation(rows, files, p),
        'sources': {os.path.basename(path): rel(path) for _m, path in files},
        'database_name': database_name,
        'database_source': database_source,
        'ext_args': threads,
    }


def sheet(name):
    return read(f'{FIX}/samplesheets/{name}').replace('${FIXTURES}', FIX)


def component_spec(metas, workflow_text, database_name, database_source, annotation_rows=None, workdir='<workDir>'):
    return {
        'workflow_text': workflow_text,
        'manifest_rows': [[n, m.get('experiment'), m.get('bioreplicate'), m['data_type'], m['group']] for n, m in metas],
        'annotation_rows': annotation_rows or [],
        'sources': {n: f'{workdir}/{n}' for n, _m in metas},
        'database_name': database_name,
        'database_source': database_source,
        'ext_args': '--threads 2',
    }


def wf_fix(n):
    return read(f'{FIX}/fragpipe_workflows/{n}')


def wf_asset(n):
    return read(f'{REPO}/assets/{n}')


TEST_FASTA = ('test.fasta', 'tests/fixtures/test.fasta')


def PHIL(src):
    return ('database_philosopher.fasta', f'<PHILOSOPHERDATABASE>/database_philosopher.fasta (from {src})')


HOSTED = ('test_decoys_contam.fasta', 'tests/fixtures/iproteomes/test_decoys_contam.fasta')
ENTRY_FASTA = ('test.fasta', 'tests/fixtures/iproteomes/test.fasta')

EB = sheet('dda_experiment_bioreplicate.csv')
BLANK = sheet('dda_experiment_blank.csv')
TMT2 = sheet('tmt_per_channel_two_plex.csv')
LC229 = wf_asset('TMT-labelcheck-229.workflow')

tests = []  # (file, test name, preferred fixture name, spec)


def add(f, name, slug, spec):
    stem = os.path.basename(f).replace('.nf.test', '')
    if stem == 'main':
        stem = {'fragpipe': 'fragpipe_wf'}[f.split('/')[2]]
    tests.append((f, name, f'{stem}__{slug}', spec))


# tests/*.nf.test (pipeline)
add('tests/fragpipe_container.nf.test', '--fragpipe_container is the container of MSFragger and IonQuant tasks', 'lfq_mbr_assets',
    pipeline_spec(EB, wf_asset('LFQ-MBR.workflow'), *TEST_FASTA))
add('tests/fragpipe_container.nf.test', 'without --fragpipe_container MSFragger and IonQuant tasks have no container', 'lfq_mbr_assets',
    pipeline_spec(EB, wf_asset('LFQ-MBR.workflow'), *TEST_FASTA))
add('tests/test_dda_lfq.nf.test', '-profile test_dda_lfq - starts with the bundled LFQ-MBR workflow and reaches IonQuant', 'lfq_mbr_philosopher_db',
    pipeline_spec(EB, wf_asset('LFQ-MBR.workflow'), *PHIL('tests/fixtures/test.fasta')))
add('tests/test_tmtq.nf.test', '-profile test_tmtq - starts with the bundled TMT6 workflow and reaches TMT-Integrator', 'tmt6_two_plex',
    pipeline_spec(sheet('tmt6_two_plex.csv'), wf_asset('TMT6.workflow'), *PHIL('tests/fixtures/test.fasta')))

# workflows/tests/fragpipe_parity.nf.test
F = 'workflows/tests/fragpipe_parity.nf.test'
add(F, 'RC-01 grouping - experiment + bioreplicate gives one filter per experiment_bioreplicate', 'lfq_mbr_experiment_bioreplicate',
    pipeline_spec(EB, wf_fix('LFQ-MBR.workflow'), *TEST_FASTA))
add(F, 'RC-04 grouping - SpecLibGen/DIA-NN workflow puts named experiments in the single group', 'dia_speclib_quant_single_group',
    pipeline_spec(EB, wf_fix('DIA_SpecLib_Quant.workflow'), *TEST_FASTA))
add(F, 'RC-01 grouping - blank experiment and bioreplicate with MBR fails like FragPipe', 'lfq_mbr_blank_experiment',
    pipeline_spec(BLANK, wf_fix('LFQ-MBR.workflow'), *TEST_FASTA))
add(F, 'RC-01 grouping - blank experiment and bioreplicate without MBR runs IonQuant on the root group', 'lfq_nombr_blank_experiment',
    pipeline_spec(BLANK, wf_fix('LFQ-noMBR.workflow'), *TEST_FASTA))
add(F, 'RC-02/RC-03 TMT per-channel samplesheet - one search per raw file, one annotation per plex', 'tmt10_per_channel_two_plex',
    pipeline_spec(TMT2, wf_fix('TMT10.workflow'), *TEST_FASTA))
add(F, 'RC-01/RC-02 TMT blank experiment - group "" through IonQuant and TMT-Integrator', 'tmt10_blank_experiment',
    pipeline_spec(sheet('tmt_blank_experiment.csv'), wf_fix('TMT10.workflow'), *TEST_FASTA))

# workflows/tests/fragpipe_run_files.nf.test
F = 'workflows/tests/fragpipe_run_files.nf.test'
add(F, 'EX-02 - pipeline mode writes fragpipe-files.fp-manifest, fragpipe.workflow and sdrf.tsv', 'lfq_mbr',
    pipeline_spec(EB, wf_fix('LFQ-MBR.workflow'), *TEST_FASTA))
add(F, 'Write sub mzML - _sub.mzML files and the second-pass workflow and manifest in fragpipe/', 'lfq_mbr_write_sub_mzml',
    pipeline_spec(EB, wf_fix('LFQ-MBR.workflow').replace('tab-run.write_sub_mzml=false', 'tab-run.write_sub_mzml=true'), *TEST_FASTA))
add(F, 'EX-02 - sdrf.tsv maps a .raw input to the file MSFragger searched', 'lfq_nombr_raw_and_mzml',
    pipeline_spec(f'sample,spectra,data_type\nthermo,{FIX}/spectra/thermo.raw,DDA\nctrl_1,{FIX}/spectra/ctrl_1.mzML,DDA\n', wf_fix('LFQ-noMBR.workflow'), *TEST_FASTA))
add(F, 'EX-02 - split-database search writes sdrf.tsv from the calibration log', 'lfq_mbr_slice_db2',
    pipeline_spec(EB, wf_fix('LFQ-MBR.workflow').replace('msfragger.misc.slice-db=1', 'msfragger.misc.slice-db=2'), *TEST_FASTA))

# workflows/tests/fragpipe_database_parity.nf.test (one test block, four iterations)
F = 'workflows/tests/fragpipe_database_parity.nf.test'
for decoy, contam in [(False, False), (True, False), (True, True), (False, True)]:
    db = PHIL('tests/fixtures/test.fasta') if (decoy or contam) else TEST_FASTA
    add(F, f'RC-07 database - decoy={str(decoy).lower()} contaminants={str(contam).lower()}',
        'lfq_mbr_database_update' if (decoy or contam) else 'lfq_mbr_database_as_given', pipeline_spec(EB, wf_fix('LFQ-MBR.workflow'), *db))

# workflows/tests/fragpipe_diapasef.nf.test: two empty .d dirs made in outputDir, DIA-typed
F = 'workflows/tests/fragpipe_diapasef.nf.test'
spec = pipeline_spec('sample,spectra,experiment,data_type\nrun_1,<outputDir>/run_1.d,exp,DIA\nrun_2,<outputDir>/run_2.d,exp,DIA\n',
                     wf_fix('DIA_SpecLib_Quant_diaPASEF.workflow'), *TEST_FASTA)
add(F, 'diaPASEF - DIA .d files go through DIATRACER, MSFragger searches the diaTracer mzML', 'diapasef_two_d', spec)

# workflows/tests/fragpipe_iproteomes.nf.test (success tests)
F = 'workflows/tests/fragpipe_iproteomes.nf.test'
add(F, 'the hosted FragPipe database is searched as given', 'lfq_mbr_hosted_db', pipeline_spec(EB, wf_fix('LFQ-MBR.workflow'), *HOSTED))
add(F, '--add_philosopher_decoy updates the plain iProteome fasta', 'lfq_mbr_iproteome_update',
    pipeline_spec(EB, wf_fix('LFQ-MBR.workflow'), *PHIL('tests/fixtures/iproteomes/test.fasta')))
add(F, 'an entry without fragpipe_fasta searches its fasta as given', 'lfq_mbr_iproteome_fasta', pipeline_spec(EB, wf_fix('LFQ-MBR.workflow'), *ENTRY_FASTA))
add(F, '--database wins over the iProteome', 'lfq_mbr_database_wins', pipeline_spec(EB, wf_fix('LFQ-MBR.workflow'), *TEST_FASTA))

# workflows/tests/tmt_labelcheck_parity.nf.test
F = 'workflows/tests/tmt_labelcheck_parity.nf.test'
add(F, 'RC-03 TMT per-channel samplesheet - one MSFRAGGER task per raw file', 'labelcheck229_two_plex', pipeline_spec(TMT2, LC229, *TEST_FASTA))
add(F, 'save-sdrf=false - the label check writes no sdrf.tsv', 'labelcheck229_no_sdrf',
    pipeline_spec(TMT2, LC229 + '\nworkflow.misc.save-sdrf=false\n', *TEST_FASTA))
add(F, 'Write sub mzML - label check writes _sub.mzML and the second-pass files to fragpipe/', 'labelcheck229_write_sub_mzml',
    pipeline_spec(TMT2, LC229 + 'tab-run.write_sub_mzml=true\n', *TEST_FASTA))

# workflows/tests/tmt_labelcheck_iproteomes.nf.test: the label check's TMT-labelcheck-229 workflow (tmt_type TMT10)
F = 'workflows/tests/tmt_labelcheck_iproteomes.nf.test'
add(F, 'the hosted FragPipe database is searched as given', 'labelcheck229_two_plex_hosted_db', pipeline_spec(TMT2, LC229, *HOSTED))
add(F, '--database wins over the iProteome', 'labelcheck229_two_plex', pipeline_spec(TMT2, LC229, *TEST_FASTA))

# subworkflows/local/fragpipe/tests/main.nf.test
F = 'subworkflows/local/fragpipe/tests/main.nf.test'
inline = re.findall(r"workflow_file\.text = '''\n(.*?)'''\.trim\(\)", read(f'{REPO}/{F}'), re.S)
assert len(inline) == 4, len(inline)
inline = [t.strip() for t in inline]
ts = ('test_sample.mzML', {'experiment': 'test_sample', 'group': 'test_sample', 'data_type': 'DDA'})
TDB = ('test_db_philosopher.fasta', '<PHILOSOPHERDATABASE>/test_db_philosopher.fasta (from <workDir>/test.fasta)')
add(F, 'fragpipe_wf - stub test with LFQ workflow', 'lfq', component_spec([ts], inline[0], *TDB))
add(F, 'fragpipe_wf - stub - export matched fragments with Report off is a no-op like FragPipe', 'matched_fragments_report_off', component_spec([ts], inline[1], *TDB))
add(F, 'fragpipe_wf - stub FPOP workflow runs FPOP Quant, coADAPTr and Metaproteomics on the final tables', 'fpop_quant_meta',
    component_spec([('fpop_a.mzML', {'experiment': 'CTRL_1', 'group': 'CTRL_1', 'data_type': 'DDA'}), ('fpop_b.mzML', {'experiment': 'OX_1', 'group': 'OX_1', 'data_type': 'DDA'})],
                   wf_fix('FPOP-quant-meta.workflow'), *TDB))
add(F, 'fragpipe_wf - stub test with TMT workflow', 'tmt16_inline',
    component_spec([ts], inline[2], *TDB, annotation_rows=[['test_sample', c, 'Bridge' if i == 0 else f'S{i}'] for i, c in enumerate(TMT16)]))
add(F, 'fragpipe_wf - stub test with prebuilt philosopher database', 'lfq_prebuilt_db',
    component_spec([ts], inline[3], 'philosopher_db.fas', '<workDir>/philosopher_db.fas'))
dia = [(f'dia_run{i}.mzML', {'group': '', 'experiment': None, 'bioreplicate': None, 'data_type': 'DIA'}) for i in (1, 2)]
add(F, 'fragpipe_wf - stub test with the bundled DIA_SpecLib_Quant workflow (SpecLibGen + DIA-NN)', 'dia_speclib_quant', component_spec(dia, wf_fix('DIA_SpecLib_Quant.workflow'), *TDB))
coadaptr = '\n'.join(line for line in wf_fix('DIA_SpecLib_Quant.workflow').splitlines() if not line.startswith('fpop.coadaptr.fpop.run-fpop-coadaptr=')) + '\nfpop.coadaptr.fpop.run-fpop-coadaptr=true\n'
add(F, "fragpipe_wf - stub DIA_SpecLib_Quant with FPOP coADAPTr edits DIA-NN's report.tsv", 'dia_speclib_quant_coadaptr', component_spec(dia, coadaptr, *TDB))
glyco = wf_fix('glyco-N-LFQ.workflow')
cancer = [('test_sample.mzML', {'group': 'cancer', 'experiment': 'cancer', 'bioreplicate': None, 'data_type': 'DDA'})]
add(F, 'fragpipe_wf - stub - glyco chain: PTM-Shepherd -> MBG -> IonQuant with AppendFile', 'glyco_n_lfq_ptms_mbg',
    component_spec(cancer, glyco.replace('ptmshepherd.run-shepherd=false', 'ptmshepherd.run-shepherd=true').replace('mbg.run-mbg=false', 'mbg.run-mbg=true')
                   .replace('saintexpress.run-saint-express=true', 'saintexpress.run-saint-express=false'), *TDB))
add(F, 'fragpipe_wf - stub - glyco chain on a FreeQuant group reads the final Philosopher report', 'glyco_n_lfq_ptms_freequant',
    component_spec(cancer, glyco.replace('ptmshepherd.run-shepherd=false', 'ptmshepherd.run-shepherd=true').replace('freequant.run-freequant=false', 'freequant.run-freequant=true')
                   .replace('ionquant.run-ionquant=true', 'ionquant.run-ionquant=false').replace('saintexpress.run-saint-express=true', 'saintexpress.run-saint-express=false'), *TDB))
ab = [('glyco_a.mzML', {'experiment': 'A', 'group': 'A', 'data_type': 'DDA'}), ('glyco_b.mzML', {'experiment': 'B', 'group': 'B', 'data_type': 'DDA'})]
add(F, 'fragpipe_wf - stub glyco O-Pair workflow: PairScans and O-Pair see one LC-MS file per run', 'glyco_o_pair', component_spec(ab, wf_fix('glyco-O-Pair.workflow'), *TDB))
add(F, 'fragpipe_wf - stub glyco-N-LFQ workflow: tsv_pepXML search without PIN, PeptideProphet per file', 'glyco_n_lfq', component_spec(ab, glyco, *TDB))

# Dedupe on the keys the recorder reads; the first test's spec (and name) wins
KEYS = ['workflow_text', 'manifest_rows', 'annotation_rows', 'database_name', 'ext_args']
os.makedirs(OUT, exist_ok=True)
by_key = {}
table = []
for f, name, fixture, spec in tests:
    key = json.dumps([spec[k] for k in KEYS])
    if key not in by_key:
        by_key[key] = fixture
        with open(f'{OUT}/{fixture}.plan_inputs.json', 'w', encoding='utf-8') as fh:
            fh.write(json.dumps(spec, indent=2, ensure_ascii=False) + '\n')
    table.append((f, name, by_key[key], spec['database_name'], [r[4] for r in spec['manifest_rows']], len(spec['annotation_rows'])))
for row in table:
    print('\t'.join(str(c) for c in row))
print('fixtures:', len(by_key), file=sys.stderr)
