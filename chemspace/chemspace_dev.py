import csv, json, argparse, copy, re, os, requests

import numpy as np
from scipy.spatial import distance
from sklearn import manifold, metrics, decomposition, preprocessing

try:
    from umap import UMAP
    from MulticoreTSNE import MulticoreTSNE
except Exception as e:
    print(e)

import igraph
import jsmin

import rdkit
from rdkit import Chem, DataStructs, Geometry
from rdkit.DataStructs import cDataStructs
from rdkit.Chem import Draw, AllChem, Scaffolds, Lipinski, Crippen, rdMolDescriptors, TemplateAlign
from rdkit.Chem.Scaffolds import MurckoScaffold

PROPS_ORDER = ["mw", "hba", "hbd", "rb", "rc", "arc", "logp", "tpsa"]
        
PROP2FNC = {
    "mw": rdMolDescriptors.CalcExactMolWt,
    "hba": Lipinski.NumHAcceptors,
    "hbd": Lipinski.NumHDonors,
    "rb": Lipinski.NumRotatableBonds,
    "rc": Lipinski.RingCount,
    "arc": Lipinski.NumAromaticRings,
    "logp": Crippen.MolLogP, 
    "tpsa": rdMolDescriptors.CalcTPSA,
}

PROP2LABEL = {
    "mw": "Molecular weight",
    "hba": "H-bond acceptors",
    "hbd": "H-bond donors",
    "rb": "Rotatable bonds",
    "rc": "Rings",
    "arc": "Aromatic rings",
    "logp": "cLogP",
    "tpsa": "TPSA"
}

FP2FNC = {
    "ecfp4": lambda rdmol: AllChem.GetMorganFingerprintAsBitVect(rdmol, radius=2, nBits=1024),
    "ecfp6": lambda rdmol: AllChem.GetMorganFingerprintAsBitVect(rdmol, radius=3, nBits=1024),
    "apfp": lambda rdmol: AllChem.GetHashedAtomPairFingerprintAsBitVect(rdmol, nBits=1024),
    "ttfp": lambda rdmol: AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(rdmol, nBits=1024),
    "maccs": lambda rdmol: AllChem.GetMACCSKeysFingerprint(rdmol),
}

AVAILABLE_METRICS = ["Tanimoto", "Dice", "Cosine", "Sokal", "Russel", "RogotGoldberg", "AllBit", "Kulczynski", "McConnaughey", "Asymmetric", "BraunBlanquet"]

DATA_KEYS = {
    "default": {},
    "compressed": {
        "features": "f",
        "object_ids": "o",
        "smiles": "s",
        "label": "l"
    }
}

class ChemSpace():

    def __init__(self, category_field = False, category_field_delimiter = False, label_field = False, compound_structure_field = False, write_structures = True, fp = "ecfp4", fingerprint_field = False, metric = "Tanimoto", pcp = False, compressed_data_format=False):
        self.category_field = category_field
        self.category_field_delimiter = category_field_delimiter
        self.label_field = label_field
        self.compound_structure_field = compound_structure_field
        self.write_structures = write_structures
        self.fp = fp
        self.fingerprint_field = fingerprint_field
        self.metric = metric
        self.pcp = pcp
        self.sdf = False
        self.KEYS = DATA_KEYS["default"] if not compressed_data_format else DATA_KEYS["compressed"]

        if self.metric not in AVAILABLE_METRICS:
            raise Exception("Metric '{}' not found in available similarity metrics: {}".format(self.metric, AVAILABLE_METRICS))

        self.index2rdmol = {}
        self.index2fpobj = {}

    def read_csv(self, filename, delimiter=",", header=False, missing_value=False, remove_columns=False):
        """Reads data from the CSV file"""
        print("Reading file: {}".format(filename))

        self.filename = filename
        with open(self.filename, "r") as input_file:
            reader = csv.reader(input_file, delimiter=delimiter)
            rows = [row for row in reader]
        
        self.read_data(rows, header, missing_value, remove_columns)

    def read_sdf(self, filename):
        """Reads data from a sdf file"""
        print("Reading file: {}".format(filename))
        self.sdf = True
        self.header = False
        self.data = []
        self.filename = filename
        self.index2rdmol = {}
        self.index2fpobj = {}
        self.index2props = {}
        self.index2category = {}
        self.index2label = {}
        self.index2id = {}

        molsupplier = Chem.SDMolSupplier(str(filename))
        not_parsed = []

        for index, m in enumerate(molsupplier):
            try:
                Chem.SanitizeMol(m)
                self.index2rdmol[index] = m
                self.index2fpobj[index] = FP2FNC[self.fp](m)
                self.index2props[index] = m.GetPropsAsDict()
                self.index2id[index] = index

            except Exception as e:
                print(e)
                not_parsed.append(index)

        self.index_order = self.index2rdmol.keys()
        self.index_order.sort()
        self.index2row = {i: [] for i in self.index_order}
        self.data = self.index2row.values()
        
        if self.label_field is not False and self.label_field in self.index2props[self.index_order[0]]:
            self.index2label = {i: self.index2props[i].get(self.label_field) for i in self.index_order}
            self.index2label = {key: val if val not in [None, ""] else self.index2id[key] for key, val in self.index2label.items()}

        if self.category_field is not False and self.category_field in self.index2props[self.index_order[0]]:
            self.index2category = {i: self.index2props[i].get(self.category_field) for i in self.index_order}

        self.__create_chemspace_format__()

    def add_compounds_from_file(self, filename, delimiter=","):
        print("Reading compounds: {}".format(filename))
        self.filename = filename

        with open(self.filename, "r") as input_file:
            reader = csv.reader(input_file, delimiter=delimiter)
            rows = [row for row in reader]

        self.add_compounds(rows)

    def add_category(self, category):
        if not "categories" in self.chemical_space:
            self.chemical_space["categories"] = []

        self.chemical_space["categories"].append(category)


    def add_compounds(self, rows):
        """Reads data in a form of list of lists (tuples)"""
        self.compounds = {r[0]: r[1] for r in rows}
        self.chemical_space["compounds"] = {}
        self.__parse_compounds__()

        for key in self.chemical_space["points"]:
            if key in self.id2rdmol:
                self.chemical_space["compounds"][key] = {"structure": self.__get_compound__(key)}

    def add_smiles(self, id2smiles):
        """Reads data in a form of list of lists (tuples)"""
        self.chemical_space["compounds"] = {}

        for key, value in self.chemical_space["points"].items():
            for i in value[self.KEYS.get("object_ids", "object_ids")]:
                if i in id2smiles:
                    self.chemical_space["compounds"][key] = {self.KEYS.get("smiles", "smiles"): id2smiles[i]}

    def read_data(self, rows, header=False, missing_value=False, remove_columns=False):
        """Reads data in a form of list of lists (tuples)"""
        self.header = header
        self.missing_value = missing_value
        data_start = 0

        self.data = rows
        self.index2id = {}
        self.index2row = {}
        self.index2compound = {}
        self.index2label = {}
        self.index2category = {}

        if self.header:
            self.header = self.data[0]
            self.data = self.data[1:]

        self.index2id = {i: row[0] for i, row in enumerate(self.data)}
        
        if self.header:
            if remove_columns is not False and len(remove_columns) > 0:
                for col in remove_columns:
                    self.__remove_field__(col)

            if self.compound_structure_field and self.compound_structure_field in self.header:
                self.index2compound = self.__extract_field__(self.compound_structure_field)
                self.__read_compounds__()

            if self.label_field and self.label_field in self.header:
                self.index2label = self.__extract_field__(self.label_field)
                self.index2label = {key: val if not self.index2label.get(key, False) else self.index2label[key] for key, val in self.index2id.items()}

            if self.category_field and self.category_field in self.header:
                self.index2category = self.__extract_field__(self.category_field)

            if self.fingerprint_field and self.fingerprint_field in self.header:
                self.index2fp = self.__extract_field__(self.fingerprint_field)
                self.index2fpobj = {}

                for index, fp in self.index2fp.items():
                    self.index2fpobj[index] = self.__get_bitvect_for_fp__(fp)

            # remove ID field
            self.header.pop(0)

        self.index2row = {i: [round(float(v), 2) if v not in ["", None, "None", self.missing_value] else None for v in row[1:]] for i, row in enumerate(self.data)}
        
        if len(self.index2rdmol) > 0:
            self.index_order = list(self.index2rdmol.keys())
            self.index_order.sort()
        else:
            self.index_order = [i for i, row in enumerate(self.data)]

        self.data = [self.index2row[i] for i in self.index_order]
        
        if self.missing_value is not False and len(self.data[0]) > 0:
            self.data, self.missing_values_indexes = self.__impute_missing_values__(self.data)

        self.__create_chemspace_format__()

    def __read_compounds__(self):
        for i, smi in self.index2compound.items():
            try:
                rdmol = Chem.MolFromSmiles(smi)

                if rdmol is not None:
                    self.index2rdmol[i] = rdmol
                    self.index2fpobj[i] = FP2FNC[self.fp](rdmol)

            except Exception as e:
                print(e)

    def __remove_field__(self, field):
        if field in self.header:
            index = self.header.index(field)

            if index is not False:
                self.header.pop(index)

                for i, row in enumerate(self.data):
                    self.data[i].pop(index)

    def __extract_field__(self, field):
        index2value = {}

        if field in self.header:
            index = self.header.index(field)

            if index is not False:
                self.header.pop(index)

                for i, row in enumerate(self.data):
                    if row[index] not in [None, "", False]:
                        index2value[i] = row[index]
                    self.data[i].pop(index)

        return index2value

    def __impute_missing_values__(self, data):
        datatype2impute = {"numeric": {"strategy":"mean", 
                                        "value": lambda x: round(float(value), 3)}, 
                           "binary": {"strategy":"most_frequent", 
                                      "value": lambda x: int(value)}
                           }

        missing_values_indexes = []
        
        for i, row in enumerate(self.data):
            missing_values_indexes.append([j for j, v in enumerate(row) if v == self.missing_value])

            for j, value in enumerate(row):
                if value == self.missing_value:
                    data[i][j] = np.nan
        imputer = preprocessing.Imputer(missing_values="NaN", strategy=datatype2impute["numeric"]["strategy"])
        #error when using median strategy - minus one dimension in imputed data... omg
        imputed_data = [list(row) for row in imputer.fit_transform(self.data)]
        imputed_data = [[datatype2impute["numeric"]["value"](value) for value in row] for row in imputed_data]
        return imputed_data, missing_values_indexes

    def __return_missing_values__(self, data, missing_values_indexes):
        for i, indexes in enumerate(missing_values_indexes):
            if indexes:
                for index in indexes:
                    data[i][index] = None
        return data

    def __create_chemspace_format__(self):
        self.chemical_space = {"points": {}}

        for index in self.index_order:
            self.chemical_space["points"][index] = {self.KEYS.get("object_ids", "object_ids"): [self.index2id[index]]}

        if len(self.index2category):
            self.__parse_categories__()

        if len(self.index2label):
            for index in self.index_order:
                self.chemical_space["points"][index][self.KEYS.get("label", "label")] = self.index2label[index]
        
        for index in self.index_order:
            self.chemical_space["points"][index][self.KEYS.get("features", "features")] = copy.copy(self.index2row[index])

        if self.header:
            current_header = self.chemical_space.get("feature_names", [])
            current_header.extend(self.header)
            self.chemical_space["feature_names"] = current_header

        if len(self.index2rdmol) and self.write_structures:
            self.chemical_space["compounds"] = {}

            for index, rdmol in self.index2rdmol.items():
                if rdmol is not None:
                    self.chemical_space["compounds"][index] = {self.KEYS.get("smiles", "smiles"): Chem.MolToSmiles(rdmol, True)}

    def __parse_categories__(self):
        category2ids = {}
        
        for index, category in self.index2category.items():
            categories = [category] if self.category_field_delimiter is False else [c.strip() for c in category.split(self.category_field_delimiter)]

            for c in categories:
                if c in category2ids:
                    category2ids[c].add(self.index2id[index])
                else:
                    category2ids[c] = {self.index2id[index]}
        
        if not "categories" in self.chemical_space:
            self.chemical_space["categories"] = []

        for c, ids in category2ids.items():
            self.chemical_space["categories"].append({self.KEYS.get("label", "label"): c, "objects": list(ids)})

    def add_paths(self, paths):
        if not self.chemical_space.get("paths", False):
            self.chemical_space["paths"] = []

        self.chemical_space["paths"].extend(paths)

    def add_paths_from_file(self):
        pass

    def add_physico_chemical_properties(self):
        print("Calculating physico-chemical properties: {} compounds".format(len(self.index2rdmol)))
        self.pcp = True
        if len(self.index2rdmol):
            count = len(self.index2rdmol)
            i = 0

            id2pcp = {}
            for index, rdmol in self.index2rdmol.items():
                if i%100 == 0 or i == count:
                    print("{}/{}".format(i, count))

                id2pcp[index] = self.__get_pcp_for_rdmol__(rdmol)
                i+=1

            empty = [None for x in PROP2LABEL]
            for i, index in enumerate(self.index_order):
                
                if id2pcp.get(index, False):
                    pcps = id2pcp[index]
                else:
                    pcps = empty
                
                self.chemical_space["points"][index][self.KEYS.get("features", "features")].extend(pcps)
                self.data[i].extend(pcps)                

            current_header = self.chemical_space.get("feature_names", [])
            current_header.extend([PROP2LABEL[prop] for prop in PROPS_ORDER])
            self.chemical_space["feature_names"] = current_header
            self.original_data = copy.deepcopy(self.data)

    def __get_pcp_for_rdmol__(self, rdmol):
        return [round(PROP2FNC[prop](rdmol), 2) for prop in PROPS_ORDER]

    def __get_compound__(self, rdmol):
        if rdmol is not None:
            Chem.Kekulize(rdmol)
            AllChem.Compute2DCoords(rdmol)
            compound = {"atoms": {}}
            atoms = [a for a in rdmol.GetAtoms()]
            bond_types = []
            for i, a in enumerate(atoms, 1):
                number = a.GetIdx()
                position = rdmol.GetConformer().GetAtomPosition(number)

                compound["atoms"][number] = {
                    "bonds": {b.GetEndAtomIdx():b.GetBondTypeAsDouble() for b in a.GetBonds() if b.GetEndAtomIdx() != number},
                    "symbol": a.GetSymbol(),
                    "charge": a.GetFormalCharge(),
                    "coordinates": [round(position.x, 3), round(position.y, 3)]
                }

                bond_types.extend(compound["atoms"][number]["bonds"].values())
        else:
            compound = None

        return compound

    def normalize_data(self, feature_range=(0,1)):
        """Normalizes data to a scale from 0 to 1."""
        print("Data normalization (scale): {}".format(feature_range))

        min_max_scaler = preprocessing.MinMaxScaler(feature_range)
        self.data = min_max_scaler.fit_transform(self.data)
        self.data = [[round(v, 3) for v in row] for row in self.data]

    def __calculate_distance_matrix__(self, similarity_threshold, index_order = None, index2order = None):
        if index_order is None:
            index_order = self.index_order

        if index2order is None:
            index2order = {x: i for i, x in enumerate(self.index_order)}
        
        print("\nCalculating distance matrix: {} compounds".format(len(index_order)))
        empty = [None for x in index_order]
        self.dist_matrix = {x:copy.copy(empty) for x in index_order}
        fps_count = len(index_order)
        
        for i, index_1 in enumerate(index_order):
            self.index2edges[index_1] = []

            if i%100 == 0 or i == fps_count:
                print("{}/{}".format(i, fps_count))
            
            for j, index_2 in enumerate(index_order[i:], i):
                sim = DataStructs.FingerprintSimilarity(self.index2fpobj[index_1], self.index2fpobj[index_2], metric=getattr(DataStructs, "{}Similarity".format(self.metric)))
                dist = 1-sim
                self.dist_matrix[index_1][j] = dist
                
                if index_1 != index_2:
                    self.dist_matrix[index_2][i] = dist

                    if sim >= similarity_threshold:
                        self.edges.append((index2order[index_1], index2order[index_2], sim))
                        self.index2edges[index_1].append([index_2, round(sim, 2)])

    def __get_edges__(self, similarity_threshold=0.7, knn=2, index_order=None, index2order=None):
        self.edges = []
        self.index2edges = {}

        if index_order is None:
            index_order = self.index_order

        if index2order is None:
            index2order = {x: i for i, x in enumerate(index_order)}
        
        count = len(index_order)
        print("\nCalculating edges [similarity threshold={}]: {} compounds".format(similarity_threshold, count))

        for i, index in enumerate(index_order):
            if (i+1)%100 == 0:
                print("{}/{}".format(i, count))

            values = [[idx, v] for idx, v in zip(index_order, self.dist_matrix[index]) if idx != index]
            values.sort(key=lambda x: x[1])

            # if 1-values[1][1] >= similarity_threshold: 
            self.index2edges[index] = []

            for v in values:
                sim = 1-v[1]
                if sim >= similarity_threshold:
                    # self.edges.append((index2order[index], index2order[v[0]]))
                    self.edges.append((index2order[index], index2order[v[0]], sim))
                    # self.index2edges[index].append([v[0]])
                    self.index2edges[index].append([v[0], round(sim, 2)])

                    if len(self.index2edges[index]) == knn:
                        break
                else:
                    break

        print("EDGES: {}".format(len(self.edges)))
        
    def __convert_fps_to_bitvects__(self, fps):
        self.index2fpobj = {}
        for i, fp in enumerate(fps):
            self.index2fpobj[i] = self.__get_bitvect_for_fp__(fp)

    def __get_bitvect_for_fp__(self, fp):
        if type(fp) is list and len(fp) == 1:
            fp = fp[0]
        bitvect = cDataStructs.ExplicitBitVect(len(fp))
        on_indexes = [i for i, b in enumerate(fp) if int(b) == 1]
        bitvect.SetBitsFromList(on_indexes)
        return bitvect

    def arrange(self, by="fps", fps=[], method="pca", similarity_threshold=0.7, add_edges=False, knn=None, add_scaffolds_category=False):
        self.dist_matrix = False
        bitvects = False
        self.edges = []
        self.index2edges = {}
        self.add_scaffolds_category = add_scaffolds_category

        if type(method) is not list:
            methods = [method]
        else:
            methods = method

        for method in methods:
            if method == "csn_scaffolds":
                scaffold_index_order, scaffold_index2order = self.__arrange_by_scaffolds__()
                self.__calculate_distance_matrix__(similarity_threshold, index_order=scaffold_index_order, index2order=scaffold_index2order)

                if knn is not None:
                    self.__get_edges__(similarity_threshold=similarity_threshold, knn=knn, index_order=scaffold_index_order, index2order=scaffold_index2order)

                g = igraph.Graph(len(self.index_order))
                print("\nCalculating Chemical Space Network...")
                feature_names = ["CSN1", "CSN2"]
                self.edges.extend(self.scaffold_edges)

                for e in self.edges:
                    g.add_edge(e[0], e[1], weight=e[2])

                layout = g.layout_fruchterman_reingold(weights="weight")
                coords = layout.coords

            elif (method in ["csn", "csn_weighted", "nn", "mds"] or by == "dm") and not self.dist_matrix:
                
                if len(fps) == 0:
                    for index in self.index_order:
                        fps.append(self.index2fpobj[index])
                    bitvects = True
                        
                elif type(fps[0][1]) in [str, int] and not bitvects:
                    self.__convert_fps_to_bitvects__(fps)
                    bitvects = True
                    
                self.__calculate_distance_matrix__(similarity_threshold)
                dist_matrix = np.matrix([np.array(self.dist_matrix[index]) for index in self.index_order])

                if knn is not None:
                    self.__get_edges__(similarity_threshold=similarity_threshold, knn=knn)

                data = dist_matrix
                by = "dm"

            elif by in ["data", "fps"]:
                if by == "fps":
                    if len(fps) == 0 and len(self.index2fpobj):
                        for index in self.index_order:
                            fps.append(self.index2fpobj[index])
                    
                    data = [[int(b) for b in fp] for fp in fps]
                else:
                    data = self.data            

            if method in ["csn", "csn_weighted"]:
                g = igraph.Graph(len(self.index_order))
                print("\nCalculating Chemical Space Network...")
                feature_names = ["CSN1", "CSN2"]                    

                print("Fruchterman-Reingold Layout calculation...")

                for e in self.edges:
                    g.add_edge(e[0], e[1], weight=e[2])
                
                if method == "csn":
                    layout = g.layout_fruchterman_reingold()
                else:
                    layout = g.layout_fruchterman_reingold(weights="weight")
                    
                coords = layout.coords

            elif method == "mds":
                print("\nCalculating MDS...")
                feature_names = ["MDS1", "MDS2"]
                # sklearn implementation
                mds = manifold.MDS(n_components=2, dissimilarity='precomputed')
                coords = mds.fit_transform(data)

                # igraph implementation
                # layout = g.layout_mds(data, 2, arpack_options=igraph.ARPACKOptions(iter=1000))
                # coords = layout.coords

            elif method == "pca":
                print("\nCalculating PCA...")
                feature_names = ["PC1", "PC2"]
                pca = decomposition.PCA(n_components=2)
                coords = pca.fit_transform(data)

            elif method == "fa":
                print("\nCalculating Factor Analysis...")
                feature_names = ["FA1", "FA2"]
                fa = decomposition.FactorAnalysis(n_components=2)
                coords = fa.fit_transform(data)

            elif method == "isomap":
                print("\nCalculating Isomap...")
                feature_names = ["Isomap1", "Isomap2"]
                isomap = manifold.Isomap(n_neighbors=200, n_components=2)
                coords = isomap.fit_transform(data)

            elif method == "multicore_tsne":
                print("\nCalculating Multicore t-SNE...")
                feature_names = ["t-SNE1", "t-SNE2"]
                tsne = MulticoreTSNE(n_components=2, metric='precomputed', n_jobs=2)
                coords = tsne.fit_transform(data)

            elif method == "tsne":
                print("\nCalculating t-SNE...")
                feature_names = ["t-SNE1", "t-SNE2"]
                tsne = manifold.TSNE(n_components=2, metric='precomputed')
                coords = tsne.fit_transform(data)

            elif method == "umap":
                print("\nCalculating UMAP...")
                feature_names = ["UMAP1", "UMAP2"]
                umap = UMAP(n_neighbors=5, min_dist=1, metric="jaccard")
                coords = umap.fit_transform(data)
                coords = [[float(x[0]), float(x[1])] for x in coords]

            # elif method == "sas":
            #     print("\nCalculating SAS...")
            #     feature_names = ["Similarity", "Activity difference"]
            #     self.chemical_space = {"points": {}, "feature_names": ["SALI"]}
            #     ai = self.header.index(self.activity_field)
            #     ids = []
            #     coords = []
                
            #     for i, index_1 in enumerate(self.index_order[:-1]):
            #         for j, index_2 in enumerate(self.index_order[i:], i):
            #             if i != j:
            #                 activity_diff = round(abs(float(self.data[i][ai]) - float(self.data[j][ai])), 2)
            #                 distance = self.dist_matrix[i][j]
            #                 distance = distance if distance > 0 else 0.01

            #                 sali = round(activity_diff/distance, 2)
            #                 coord = [round(1 - self.dist_matrix[i][j], 2), activity_diff]
            #                 self.chemical_space["points"]["{}_{}".format(index_1, index_2)] = {"features": [sali]}
            #                 ids.append("{}_{}".format(index_1, index_2))
            #                 coords.append(coord)

            #     self.index_order = ids

            if method in ["csn", "csn_weighted", "csn_scaffolds", "nn"] or add_edges:
                if self.dist_matrix is False and self.edges in [False, []]:
                    self.__calculate_distance_matrix__()

                if self.edges in [False, []]:
                    if knn is None:
                        knn = len(self.index_order)
                    self.__get_edges__(similarity_threshold=similarity_threshold, knn=knn)

                
                for cid, es in self.index2edges.items():
                    if not self.chemical_space["points"][cid].get("links", False):
                        self.chemical_space["points"][cid]["links"] = []

                    for e in es:
                        self.chemical_space["points"][cid]["links"].append(e)

            index2coords = {index:coords[i] for i, index in enumerate(self.index_order)}

            for index, values in self.chemical_space["points"].items():
                if index in index2coords:
                    point_features = self.chemical_space["points"][index][self.KEYS.get("features", "features")]
                    features = [round(index2coords[index][0], 3), round(index2coords[index][1], 3)]
                    features.extend(point_features)
                    self.chemical_space["points"][index][self.KEYS.get("features", "features")] = features

                else:
                    self.chemical_space["points"].pop(index, None)

            feature_names.extend(self.chemical_space.get("feature_names", []))
            self.chemical_space["feature_names"] = feature_names

    def __arrange_by_scaffolds__(self, align_by_scaffold=True):
        self.scaffold2indexes = {}
        self.scaffold2index_orders = {}
        self.scaffold2rdmol = {}
        self.index2scaffold = {}

        for i, index in enumerate(self.index_order):
            rdmol = self.index2rdmol[index]
            scaffold = Scaffolds.MurckoScaffold.GetScaffoldForMol(rdmol)
            scaffold_smiles = Chem.MolToSmiles(scaffold, isomericSmiles=False)
            self.scaffold2rdmol[scaffold_smiles] = scaffold

            if scaffold_smiles in self.scaffold2indexes:
                self.scaffold2indexes[scaffold_smiles].append(index)
                self.scaffold2index_orders[scaffold_smiles].append(i)
            else:
                self.scaffold2indexes[scaffold_smiles] = [index]
                self.scaffold2index_orders[scaffold_smiles] = [i]

        index2order = {}
        index_order = []
        order = len(self.index_order)

        for index, scaffold in enumerate(self.scaffold2indexes.keys(), self.index_order[-1]+1):
            self.index2fpobj[index] = FP2FNC[self.fp](self.scaffold2rdmol[scaffold])
            self.index2scaffold[index] = scaffold
            self.index2rdmol[index] = self.scaffold2rdmol[scaffold]
            self.index_order.append(index)
            index_order.append(index)
            index2order[index] = order
            order+=1

        self.scaffold_edges = []

        for index_1, scaffold in self.index2scaffold.items():
            for index_2 in self.scaffold2index_orders[scaffold]:
                self.scaffold_edges.append((index2order[index_1], index_2, 1))

        self.__add_scaffolds_to_chemical_space__()

        return index_order, index2order

    def __add_scaffolds_to_chemical_space__(self):
        for index, scaffold in self.index2scaffold.items():
            print(self.scaffold2indexes[scaffold])
            self.chemical_space["points"][index] = {
                self.KEYS.get("object_ids", "object_ids"): [index],
                "links": [[x, 1] for x in self.scaffold2indexes[scaffold]],
                self.KEYS.get("label", "label"): "Scaffold {}".format(index-len(self.data)+1),
                self.KEYS.get("features", "features"): []
            }
            
            for i, f in enumerate(self.chemical_space["feature_names"]):
                values = [self.chemical_space["points"][x][self.KEYS.get("features", "features")][i] for x in self.scaffold2indexes[scaffold]]
                values = [v for v in values if v is not None]
                value = round(np.mean(values), 2) if len(values) else None
                self.chemical_space["points"][index][self.KEYS.get("features", "features")].append(value)

            self.chemical_space["compounds"][index] = {self.KEYS.get("smiles", "smiles"): scaffold, "color": "red"}

        if self.add_scaffolds_category:
            if not self.chemical_space.get("categories", False):
                self.chemical_space["categories"] = []

            self.chemical_space["categories"].append({"label": "scaffolds", "color": "gray", "points": self.index2scaffold.keys(), "shape": "circle"})

    def export_chemical_space_as_html(self, htmldir=".", ):
        """Export a simple HTML page with embedded chemical space and dependencies into a given directory."""
        if not os.path.exists(htmldir):
            os.makedirs(htmldir)

        chemspace_json = self.export_chemical_space_as_json(minify=True, dump=True)
        
        libs = [
            ("chemspace-0.2.0.min.js", "https://openscreen.cz/software/chemspace/static/js/chemspace-0.2.0.min.js"),
            ("jquery-3.3.1.min.js", "https://code.jquery.com/jquery-3.3.1.min.js"),
            ("konva.min.js", "https://cdn.rawgit.com/konvajs/konva/1.7.6/konva.min.js")
        ]
        
        js_html = []
        for l in libs:
            js_html.append("<script src='{}'></script>".format(l[0]))

        settings = {
            "target": "chemspace"
        }

        template = """<html>
        <head>
            {}
            <script>
                $(document).ready(function() {{
                    var data = {};
                    var chemspace = new ChemSpace({});
                    chemspace.read_data(data);
                    chemspace.draw();
                }});
            </script>
        </head>

        <body>
            <div id="chemspace"></div>
        </body>
        </html>""".format('\n'.join(js_html), chemspace_json, json.dumps(settings))

        
        for l in libs:
            lib, url = l
            try:
                source = requests.get(url)
                source_html = source.read()

                with open(os.path.join(htmldir, lib), "w") as output:
                    output.write(source_html)
            except requests.exceptions.RequestException as e:
                raise Exception("""
                        \nCan't download file {}.\nPlease check your internet connection and try again.\nIf the error persists there can be something wrong with the InCHlib server.\n""".format(url)
                    )

        with open(os.path.join(htmldir, "chemspace.html"), "w") as output:
            output.write(template)

    def export_chemical_space_as_json(self, filename=None, minify=False, dump=True):
        """Returns space in a JSON format or exports it to the file specified by the filename parameter."""
        space_json = self.chemical_space

        if minify:
            space_json = json.dumps(space_json)
            space_json = self.__minify_data(space_json)
        elif dump:
            space_json = json.dumps(space_json, indent=4)

        if filename:
            output = open(filename, "w")
            output.write(space_json)
        
        return space_json

    def __minify_data(self, data):
        return jsmin.jsmin(data)

def _process_(arguments):
    s = ChemSpace(
        write_structures = False if arguments.dont_write_structures else True,
        fp = arguments.fingerprint,
        category_field = arguments.category_field,
        category_field_delimiter = arguments.category_field_delimiter,
        label_field = arguments.label_field,
        compound_structure_field = arguments.compound_structure_field,
        fingerprint_field = arguments.fingerprint_field,
        metric = arguments.similarity_metric,
        compressed_data_format=arguments.compressed_data_format
    )

    if arguments.data_file.split(".")[-1].lower() == "sdf":
        s.read_sdf(arguments.data_file)
    else:
        s.read_csv(arguments.data_file, arguments.data_delimiter, arguments.data_header, arguments.missing_values, arguments.remove_columns)

    if s.compound_structure_field is not False or s.sdf == True:
        if arguments.physico_chemical_properties:
            s.add_physico_chemical_properties()

    if arguments.normalize:
        s.normalize_data()

    if arguments.arrange_by:
        s.arrange(
            method=arguments.dimensional_reduction_method,
            similarity_threshold=float(arguments.compound_similarity_threshold),
            add_edges=arguments.add_edges,
            by=arguments.arrange_by,
            knn=arguments.knn,
            add_scaffolds_category=arguments.add_scaffolds_category
        )
    
    if arguments.html_dir:
        s.export_chemical_space_as_html(arguments.html_dir)
    elif arguments.output_file:
        s.export_chemical_space_as_json(arguments.output_file, minify=arguments.minify_output)
    else:
        print(s.export_chemical_space_as_json(minify=arguments.minify_output))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("data_file", type=str, help="csv(text) data file with delimited values or a sdf file")
    parser.add_argument("-dh", "--data_header", default=False, help="whether the first row of data file is a header", action="store_true")
    parser.add_argument("-dd", "--data_delimiter", type=str, default=",", help="delimiter of values in data file")
    parser.add_argument("-o", "--output_file", type=str, help="the name of output file")
    parser.add_argument("-fpf", "--fingerprint_field", type=str, default=False, help="set a fingerprint field name in case it is in the data file")
    parser.add_argument("-cf", "--category_field", type=str, default=False, help="set a category field name in case it is in the data file")
    parser.add_argument("-cfd", "--category_field_delimiter", type=str, default=False, help="a category field delimiter")
    parser.add_argument("-lf", "--label_field", type=str, default=False, help="set a label field name in case it is in the data file")
    parser.add_argument("-af", "--activity_field", type=str, default=False, help="set an activity field name in case it is in the data file")
    parser.add_argument("-csf", "--compound_structure_field", type=str, default=False, help="the name of a column with a compound structure")
    parser.add_argument("-fp", "--fingerprint", type=str, default="ecfp4", help="fingerprint used for a compound representation (ecfp4, ecfp6, maccs, topological, atom_pairs)")
    # parser.add_argument("-c", "--compounds", type=str, default=None, help="csv(text) compound file with delimited values in a form: id,smiles")
    # parser.add_argument("-cd", "--compounds_delimiter", type=str, default=",", help="delimiter of values in compound file")
    parser.add_argument("-arr", "--arrange_by", default="data", help="arrange data by compound structures (distance matrix) or by input data (data/fps)", type=str)
    parser.add_argument("-asc", "--add_scaffolds_category", default=False, help="add scaffolds as a category (only for arrange by csn_scaffolds)", type=str)
    parser.add_argument("-cst", "--compound_similarity_threshold", default=0.7, help="compound similarity threshold")
    parser.add_argument("-drm", "--dimensional_reduction_method", nargs='+', type=str, default="pca", help="which method use for dimensional reduction (pca/isomap/csn)")
    parser.add_argument("-dws", "--dont_write_structures", default=False, help="dont write structures to output file", action="store_true")
    parser.add_argument("-min", "--minify_output", default=False, help="minify the JSON output format", action="store_true")
    parser.add_argument("-html", "--html_dir", type=str, default=False, help="the directory to store HTML page with dependencies")
    parser.add_argument("-pcp", "--physico_chemical_properties", default=False, help="calculate basic phyisico-chemical properties and add them ass features", action='store_true')
    parser.add_argument("-edges", "--add_edges", default=False, help="add edges based on compound similarity to the graph", action='store_true')
    parser.add_argument("-n", "--normalize", default=False, help="normalize data to [0, 1] range", action="store_true")
    parser.add_argument("-mv", "--missing_values", type=str, default=False, help="define the string representating missing values in the data")
    parser.add_argument("-knn", "--knn", type=int, default=None, help="the number of neighbours (k) used for the construction of csn using the nn method")
    parser.add_argument("-sm", "--similarity_metric", type=str, default="Tanimoto", help="similarity metric")
    parser.add_argument('-rmc','--remove_columns', nargs='+', default=False, help='columns in data that should not be used')
    parser.add_argument('-cdf','--compressed_data_format', nargs='+', default=False, help='use shorter data keys')
    
    args = parser.parse_args()
    _process_(args)
