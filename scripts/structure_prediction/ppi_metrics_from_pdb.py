import numpy as np
import sys
from collections import defaultdict

def parse_atm_record(line):
    return {
        'atom_name': line[12:16].strip(),
        'res_no': int(line[22:26]),
        'chain': line[21],
        'x': float(line[30:38]),
        'y': float(line[38:46]),
        'z': float(line[46:54]),
        'b_factor': float(line[60:66])
    }

def extract_cb_and_plddt(pdb_file):
    cb_coords_A, plddt_A = [], []
    cb_coords_B, plddt_B = [], []
    seen_cb = set()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                record = parse_atm_record(line)
                if record['atom_name'] != 'CB':
                    continue
                key = (record['chain'], record['res_no'])
                if key in seen_cb:  # Avoid duplicate atoms
                    continue
                seen_cb.add(key)

                coord = np.array([record['x'], record['y'], record['z']])
                if record['chain'] == 'A':
                    cb_coords_A.append(coord)
                    plddt_A.append(record['b_factor'])
                elif record['chain'] == 'B':
                    cb_coords_B.append(coord)
                    plddt_B.append(record['b_factor'])
    return np.array(cb_coords_A), np.array(plddt_A), np.array(cb_coords_B), np.array(plddt_B)

def compute_interface_contacts(cb_A, cb_B, threshold=8.0):
    if cb_A.size == 0 or cb_B.size == 0:
        return np.array([])
    dists = np.linalg.norm(cb_A[:, None, :] - cb_B[None, :, :], axis=-1)
    contacts = np.argwhere(dists <= threshold)
    return contacts

def compute_pdockq(avg_if_plddt, n_contacts):
    if n_contacts < 1:
        return 0
    x = avg_if_plddt * np.log10(n_contacts)
    return 0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python ppi_metrics_from_pdb.py <input.pdb>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    cb_A, plddt_A, cb_B, plddt_B = extract_cb_and_plddt(pdb_file)
    contacts = compute_interface_contacts(cb_A, cb_B)

    if contacts.shape[0] < 1:
        print("No interface contacts found.")
        print("num_contacts: 0")
        print("avg_if_plddt: 0.0")
        print("pDockQ: 0.0")
    else:
        # Interface residues
        res_A = np.unique(contacts[:, 0])
        res_B = np.unique(contacts[:, 1])

        avg_if_plddt = np.mean(np.concatenate([plddt_A[res_A], plddt_B[res_B]]))
        n_contacts = contacts.shape[0]
        pdockq = compute_pdockq(avg_if_plddt, n_contacts)
        mod_name= pdb_file.replace('_model.pdb','')
        # print(f"num_contacts: {n_contacts}")
        # print(f"avg_if_plddt: {avg_if_plddt:.2f}")
        # print(f"pDockQ: {pdockq:.3f}")
        print(f"{mod_name}\t{n_contacts}\t{avg_if_plddt:.2f}\t{pdockq:.3f}")

