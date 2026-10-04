"""Check that the single-zone Nagao test model reproduces grid model 870
(log phi(H) = 19, log n_H = 10, Z = 1) of the paper's Cloudy grid.

Compares (1) every line flux in the saved line list and (2) the gas-phase
chemical composition block printed in the .out files.
"""
import numpy as np

GRID_DIR = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/'
TEST_DIR = '/Users/jiamuh/c23.01/my_models/sed_test_ngc5548/'
GRID_STEM = 'strong_LOC_varym_N25_v100_lineflux'
TEST_STEM = 'nagao_n10_phi19'
GRID_INDEX = 870

KEY_LINES = ['H  1 1215.67A', 'blnd 1549.00A', 'blnd 1909.00A',
             'blnd 2798.00A', 'H  1 4861.32A', 'H  1 6562.80A']


def read_test_linelist():
    with open(TEST_DIR + TEST_STEM + '_LineList_BLR_Fe2_flux.txt') as f:
        header = f.readline().rstrip('\n').split('\t')[1:]
        values = f.readline().rstrip('\n').split('\t')[1:]
    return header, np.array(values, dtype=float)


def read_grid_linelist(index):
    delimiter = 'GRID_DELIMIT -- grid%09d' % index
    with open(GRID_DIR + GRID_STEM + '_LineList_BLR_Fe2_flux.txt') as f:
        lines = f.read().split('\n')
    header = lines[0].split('\t')[1:]
    row = next(i for i, l in enumerate(lines) if l.startswith('#') and delimiter in l)
    values = lines[row - 1].split('\t')[1:]
    return header, np.array(values, dtype=float)


def composition_blocks(path):
    with open(path) as f:
        lines = f.readlines()
    return [lines[i + 1] + lines[i + 2] for i, l in enumerate(lines)
            if 'Gas Phase Chemical Composition' in l]


test_header, test_flux = read_test_linelist()
grid_header, grid_flux = read_grid_linelist(GRID_INDEX)
print('Line list headers identical:', test_header == grid_header)

both_positive = (test_flux > 0) & (grid_flux > 0)
one_zero = (test_flux > 0) != (grid_flux > 0)
log_diff = np.log10(test_flux[both_positive] / grid_flux[both_positive])
print('Lines with positive flux in both runs: %d' % both_positive.sum())
print('Lines positive in only one run: %d' % one_zero.sum())
worst = np.flatnonzero(both_positive)[np.argmax(np.abs(log_diff))]
print('Largest |log10(test / grid)| over all lines: %.2e dex (%s, test %.3e, grid %.3e)'
      % (np.abs(log_diff).max(), test_header[worst], test_flux[worst], grid_flux[worst]))


def column(name):
    return next(i for i, h in enumerate(test_header) if h.strip().lower() == name.lower())


print('\n%-16s %14s %14s %14s' % ('line', 'test flux', 'grid flux', 'log ratio'))
for name in KEY_LINES:
    i = column(name)
    print('%-16s %14.4e %14.4e %14.2e' % (name, test_flux[i], grid_flux[i],
                                          np.log10(test_flux[i] / grid_flux[i])))

hbeta = column('H  1 4861.32A')
print('\n%-16s %14s %14s %14s' % ('line / H-beta', 'test ratio', 'grid ratio', 'log difference'))
for name in KEY_LINES:
    i = column(name)
    test_ratio = test_flux[i] / test_flux[hbeta]
    grid_ratio = grid_flux[i] / grid_flux[hbeta]
    print('%-16s %14.4f %14.4f %14.2e' % (name, test_ratio, grid_ratio,
                                          np.log10(test_ratio / grid_ratio)))

test_comp = composition_blocks(TEST_DIR + TEST_STEM + '.out')
grid_comp = composition_blocks(GRID_DIR + GRID_STEM + '.out')
print('\nComposition blocks found: test %d, grid %d' % (len(test_comp), len(grid_comp)))
print('Test composition:\n' + test_comp[0])
print('Grid model %d composition:\n' % GRID_INDEX + grid_comp[GRID_INDEX])
print('Composition identical:', test_comp[0] == grid_comp[GRID_INDEX])
