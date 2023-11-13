import module_io.read_wfc as read_wfc
import matplotlib.pyplot as plt

def read_wfc_test():
    psi = read_wfc.read_lowf('./module_io/test/support', 8, 14, 26)

    ncols = 4
    nrows = psi.get_nkpts()//ncols
    if psi.get_nkpts()%ncols != 0:
        nrows += 1

    overlap = []
    for ikpt in range(psi.get_nkpts()):
        overlap.append(psi.overlap_with(ikpt, psi))

    fig, axs = plt.subplots(nrows, ncols, figsize=(20, 20))
    for ikpt in range(psi.get_nkpts()):
        irow = ikpt//ncols
        icol = ikpt%ncols
        axs[irow, icol].imshow(overlap[ikpt].real, vmin=-1, vmax=1)
        axs[irow, icol].set_title(f'k = {ikpt+1}')

    plt.colorbar(axs[0, 0].imshow(overlap[0].real, vmin=-1, vmax=1), ax=axs, orientation='horizontal', fraction=.1)
    plt.show()
    