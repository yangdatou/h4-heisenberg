import numpy, scipy, sys
from pyscf import gto, scf
from pyscf.scf import stability

from pyscf.tools.dump_mat import dump_rec

s1 = numpy.array([[0.0,   1.0 ], [1.0,  0.0]])
s2 = numpy.array([[0.0,  -1.0j], [1.0j, 0.0]])
s3 = numpy.array([[1.0,   0.0 ], [0.0, -1.0]])
pauli_matrices = numpy.asarray([s1, s2, s3]).reshape(3, 2, 2)

r = 1.6
mol = gto.Mole()
mol.atom = f'''
    H 0 0 { 2.0*r: 6.4f}; H 0 0 { 1.0*r: 6.4f};
    H 0 0 {-2.0*r: 6.4f}; H 0 0 {-1.0*r: 6.4f};
'''
mol.spin = 0

mol.basis = 'sto-3g'
mol.build()

dm_init_1 = (numpy.diag([1, 0, 0, 1]), numpy.diag([0, 1, 1, 0]))
dm_init_2 = (numpy.diag([1, 0, 1, 0]), numpy.diag([0, 1, 0, 1]))

def build_ghf(mol, dm_init=None):
    nao = mol.nao_nr()
    ao_labels = mol.ao_labels()

    atm_ao_list = []
    for iatm in range(mol.natm):
        iatm_ao_list = []
        for mu, ao_label_mu in enumerate(ao_labels):
            ao_label_mu_split_0 = int(ao_label_mu.split()[0])
            ao_label_mu_split_1 = ao_label_mu.split()[1]
            
            if ao_label_mu_split_0 == iatm:
                assert ao_label_mu_split_1 == mol.atom_pure_symbol(iatm)
                iatm_ao_list.append(mu)
            
        atm_ao_list.append(iatm_ao_list)

    uhf = scf.UHF(mol)
    uhf.verbose = 4
    ovlp = uhf.get_ovlp()
    assert ovlp.shape == (nao, nao)

    e  = uhf.kernel(dm0=dm_init)
    dm = uhf.make_rdm1()
    assert uhf.converged

    ghf = scf.GHF(mol)
    ghf.convert_from_(uhf)

    dm = ghf.make_rdm1()
    dm_aa = dm[:nao,:nao]
    dm_bb = dm[nao:,nao:]
    dm_ab = dm[:nao,nao:]
    dm_ba = dm[nao:,:nao]

    dump_rec(sys.stdout, dm)

    # GHF density matrix in two component spinor basis
    dm_s = [[dm_aa, dm_ab], [dm_ba, dm_bb]]
    dm_s = numpy.asarray(dm_s) # .reshape(2, 2, nao, nao)
    dm_x = numpy.einsum("abmn,xab->xmn", dm_s, pauli_matrices)

    spin_on_atoms = []

    for iatm in range(mol.natm):
        iatm_ao_list      = atm_ao_list[iatm]
        ao_idx_on_atom    = numpy.ix_(iatm_ao_list, iatm_ao_list)
        spin_on_atom_x    = []

        pop_on_atom  = numpy.einsum("mn,nm->", dm_aa[iatm_ao_list, :], ovlp[:, iatm_ao_list])
        pop_on_atom -= numpy.einsum("mn,nm->", dm_bb[iatm_ao_list, :], ovlp[:, iatm_ao_list])
        print_info = f"{iatm: 3d} afm = {pop_on_atom: 6.4f}, spin = ["

        for x in range(3):
            spin_on_atom = numpy.einsum("mn,nm->", dm_x[x][iatm_ao_list, :], ovlp[:, iatm_ao_list])
            spin_on_atom_x.append(spin_on_atom)

            print_info += f" ({spin_on_atom.real: 6.4f}, {spin_on_atom.imag: 6.4f})"

            if x != 2:
                print_info += ", "
            else:
                print_info += " ]"

        spin_on_atom_x = numpy.asarray(spin_on_atom_x)
        print(f"spin_on_atom_x = {print_info}, norm = {numpy.linalg.norm(spin_on_atom_x.real): 6.4f}")


build_ghf(mol, dm_init_1)
build_ghf(mol, dm_init_2)

