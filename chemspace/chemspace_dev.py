import csv, argparse, copy, re, os, requests, time
import simplejson as json
from collections import defaultdict

import numpy as np
from scipy.spatial import distance
from sklearn import manifold, metrics, decomposition, preprocessing
from sklearn.impute import SimpleImputer

import tmap

try:
    from MulticoreTSNE import MulticoreTSNE
except Exception as e:
    print(e)

try:
    from umap import UMAP
except Exception as e:
    print(e)

import igraph
import jsmin

import rdkit
from rdkit import Chem, DataStructs, Geometry
from rdkit.DataStructs import cDataStructs
from rdkit.Chem import Draw, AllChem, Scaffolds, Lipinski,\
    Crippen, rdMolDescriptors, TemplateAlign, QED
from rdkit.Chem.Scaffolds import MurckoScaffold

PROPS_ORDER = ["mw", "hba", "hbd", "rb", "rc", "arc", "logp", "tpsa", "fcsp3", "ncc", "qed"]
        
PROP2FNC = {
    "mw": rdMolDescriptors.CalcExactMolWt,
    "hba": Lipinski.NumHAcceptors,
    "hbd": Lipinski.NumHDonors,
    "rb": Lipinski.NumRotatableBonds,
    "rc": Lipinski.RingCount,
    "arc": Lipinski.NumAromaticRings,
    "logp": Crippen.MolLogP, 
    "tpsa": rdMolDescriptors.CalcTPSA,
    "fcsp3": rdMolDescriptors.CalcFractionCSP3,
    "ncc": lambda m: len(Chem.FindMolChiralCenters(m, includeUnassigned=True)),
    "qed": QED.default
}

PROP2LABEL = {
    "mw": "Molecular weight",
    "hba": "H-bond acceptors",
    "hbd": "H-bond donors",
    "rb": "Rotatable bonds",
    "rc": "Rings",
    "arc": "Aromatic rings",
    "logp": "cLogP",
    "tpsa": "TPSA",
    "fcsp3": "Fraction of csp3",
    "ncc": "Chiral centers",
    "qed": "QED",
}

FP2FNC = {
    "ecfp4": lambda rdmol: AllChem.GetMorganFingerprintAsBitVect(rdmol, radius=2, nBits=1024),
    "ecfp6": lambda rdmol: AllChem.GetMorganFingerprintAsBitVect(rdmol, radius=3, nBits=1024),
    "atompairs": lambda rdmol: AllChem.GetHashedAtomPairFingerprintAsBitVect(rdmol, nBits=1024),
    "torsion": lambda rdmol: AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(rdmol, nBits=1024),
    "maccs": lambda rdmol: AllChem.GetMACCSKeysFingerprint(rdmol),
}

METHODS = {
    "pca": {
        "dm": False,
        "edges": False,
        "label": "Principal Components Analysis",
        "dim_label": "PCA"
    },
    "mds": {
        "dm": False,
        "edges": False,
        "label": "Multi-dimensional Scaling",
        "dim_label": "MDS"
    },
    "umap": {
        "dm": False,
        "edges": False,
        "label": "UMAP",
        "dim_label": "UMAP"
    },
    "tsne": {
        "dm": False,
        "edges": False,
        "label": "t-SNE",
        "dim_label": "t-SNE"
    },
    "multicore_tsne": {
        "dm": False,
        "edges": False,
        "label": "t-SNE",
        "dim_label": "t-SNE"
    },
    "csn": {
        "dm": True,
        "edges": True,
        "label": "Chemical Space Network",
        "dim_label": "CNS"
    },
    "csn_scaffolds": {
        "dm": True,
        "edges": True,
        "label": "Scaffold Chemical Space Network",
        "dim_label": "Scaffold CSN"
    },
    "nn": {
        "dm": True,
        "edges": False,
        "label": "Nearest Neighbours",
        "dim_label": "NN"
    },
    "mst": {
        "dm": True,
        "edges": True,
        "label": "Minimum Spanning Tree",
        "dim_label": "MST"
    },
    "mst_scaffolds": {
        "dm": True,
        "edges": True,
        "label": "Minimum Scaffold Spanning Tree",
        "dim_label": "Scaffold MST"
    },
}

AVAILABLE_METRICS = [
    "Tanimoto", "Dice", "Cosine", "Sokal", "Russel", "RogotGoldberg", 
    "AllBit", "Kulczynski", "McConnaughey", "Asymmetric", "BraunBlanquet"
]

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

    def __init__(
            self, category_field = False, category_field_delimiter = False, label_field = False, 
            compound_structure_field = False, write_structures = True, fp = "ecfp4", 
            fingerprint_field = False, metric = "Tanimoto", pcp = False, compressed_data_format=False, 
            round_values=False, data_as_features=True, n_jobs=1, keep_unparsable_structures=False
        ):
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
        self.round_values = round_values
        self.data_as_features = data_as_features
        self.keep_unparsable_structures = keep_unparsable_structures
        self.n_jobs = n_jobs
        self.KEYS = DATA_KEYS["default"] if not compressed_data_format else DATA_KEYS["compressed"]

        # if self.metric not in AVAILABLE_METRICS:
        #     raise Exception("Metric '{}' not found in available similarity metrics: {}".format(self.metric, AVAILABLE_METRICS))

        self.index2rdmol = {}
        self.index2fpobj = {}

    def read_file(self, filename, delimiter=",", header=False, missing_value=False, remove_columns=False):
        ext = filename.split(".")[-1].lower()

        if ext in ["csv", "sdf"]:
            if ext == "csv":
                self.read_csv(filename, delimiter, header, missing_value, remove_columns)
            else:
                self.read_sdf(filename)
        else:
            raise Exception(f"Input file must be of type: csv, sdf. Input file: {filename}")

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
        mi = 0

        for index, m in enumerate(molsupplier):
            try:
                Chem.SanitizeMol(m)
                self.index2rdmol[mi] = m
                self.index2fpobj[mi] = FP2FNC[self.fp](m)
                self.index2props[mi] = m.GetPropsAsDict()
                self.index2id[mi] = mi
                mi += 1

            except Exception as e:
                print(e)
                not_parsed.append(index)

        self.index_order = list(self.index2rdmol.keys())
        self.index_order.sort()
        self.index2row = {i: [] for i in self.index_order}
        self.data = list(self.index2row.values())
        
        if self.label_field is not False and self.label_field in self.index2props[self.index_order[0]]:
            self.index2label = {i: self.index2props[i].get(self.label_field) for i in self.index_order}
            self.index2label = {key: val if val not in [None, ""] else self.index2id[key] for key, val in self.index2label.items()}

        if self.category_field is not False and self.category_field in self.index2props[self.index_order[0]]:
            self.index2category = {i: self.index2props[i].get(self.category_field) for i in self.index_order}

        self._create_chemspace_format()

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


    # def add_compounds(self, rows):
    #     """Reads data in a form of list of lists (tuples)"""
    #     self.compounds = {r[0]: r[1] for r in rows}
    #     self.chemical_space["compounds"] = {}
    #     self._parse_compounds_()

    #     for key in self.chemical_space["points"]:
    #         if key in self.id2rdmol:
    #             self.chemical_space["compounds"][key] = {"structure": self._get_compound_(key)}

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
                    self._remove_field(col)

            if self.compound_structure_field and self.compound_structure_field in self.header:
                print("Extracting structure field")
                self.index2compound = self._extract_field(self.compound_structure_field)
                self._read_compounds()

            if self.label_field and self.label_field in self.header:
                print("Extracting label field")
                self.index2label = self._extract_field(self.label_field)
                self.index2label = {key: val if not self.index2label.get(key, False) else self.index2label[key] for key, val in self.index2id.items()}

            if self.category_field and self.category_field in self.header:
                print("Extracting category field")
                self.index2category = self._extract_field(self.category_field)

            if self.fingerprint_field and self.fingerprint_field in self.header:
                self.index2fp = self._extract_field(self.fingerprint_field)
                self.index2fpobj = {}

                for index, fp in self.index2fp.items():
                    self.index2fpobj[index] = self._get_bitvect_for_fp(fp)

            # remove ID field
            self.header.pop(0)

        if self.round_values is not False:
            self.index2row = {i: [round(float(v), self.round_values) if v not in ["", None, "None", self.missing_value] else None for v in row[1:]] for i, row in enumerate(self.data)}
        else:
            self.index2row = {i: [float(v) if v not in ["", None, "None", self.missing_value] else None for v in row[1:]] for i, row in enumerate(self.data)}
        
        if len(self.index2rdmol) > 0 and not self.keep_unparsable_structures:
            self.index_order = list(self.index2rdmol.keys())
            self.index_order.sort()
        else:
            self.index_order = [i for i, row in enumerate(self.data)]
        
        self.original_data = copy.deepcopy(list(self.index2row.values()))
        self.data = [self.index2row[i] for i in self.index_order]
        
        if self.missing_value is not False and len(self.data[0]) > 0:
            self.data, self.missing_values_indexes = self._impute_missing_values(self.data)
        
        self._create_chemspace_format()

    def _read_compounds(self):
        empty = Chem.MolFromSmiles("")
        empty_fp = FP2FNC[self.fp](empty)

        for i, smi in self.index2compound.items():
            try:
                rdmol = Chem.MolFromSmiles(smi)

                if rdmol is not None:
                    self.index2rdmol[i] = rdmol
                    self.index2fpobj[i] = FP2FNC[self.fp](rdmol)

            except Exception as e:
                self.index2rdmol[i] = empty
                self.index2fpobj[i] = empty_fp
                print(e)

    def _remove_field(self, field):
        if field in self.header:
            index = self.header.index(field)

            if index is not False:
                self.header.pop(index)

                for i, row in enumerate(self.data):
                    self.data[i].pop(index)

    def _extract_field(self, field):
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

    def _impute_missing_values(self, data):
        datatype2impute = {
            "numeric": {
                "strategy":"mean", 
                "value": lambda x: round(float(value), 3)
            }, 
            "binary": {
                "strategy":"most_frequent", 
                "value": lambda x: int(value)
            }
        }

        missing_values_indexes = []
        
        for i, row in enumerate(self.data):
            missing_values_indexes.append([j for j, v in enumerate(row) if v == self.missing_value])

            for j, value in enumerate(row):
                if value == self.missing_value:
                    data[i][j] = np.nan

        imputer = SimpleImputer(missing_values=np.nan, strategy=datatype2impute["numeric"]["strategy"], keep_empty_features=True)
        #error when using median strategy - minus one dimension in imputed data... omg
        imputed_data = [list(row) for row in imputer.fit_transform(self.data)]
        # imputed_data = [[datatype2impute["numeric"]["value"](value) for value in row] for row in imputed_data]
        return imputed_data, missing_values_indexes

    def _return_missing_values(self, data, missing_values_indexes):
        for i, indexes in enumerate(missing_values_indexes):
            if indexes:
                for index in indexes:
                    data[i][index] = None

        return data

    def _create_chemspace_format(self):
        self.chemical_space = {"points": {}}

        for index in self.index_order:
            self.chemical_space["points"][index] = {self.KEYS.get("object_ids", "object_ids"): [str(self.index2id[index])]}

        if len(self.index2category):
            self._parse_categories()

        if len(self.index2label):
            for index in self.index_order:
                self.chemical_space["points"][index][self.KEYS.get("label", "label")] = self.index2label[index]
        
        if self.data_as_features:
            for index in self.index_order:
                self.chemical_space["points"][index][self.KEYS.get("features", "features")] = copy.copy(self.index2row[index])

        else:
            for index in self.index_order:
                self.chemical_space["points"][index][self.KEYS.get("features", "features")] = []

        if self.header and self.data_as_features:
            current_header = self.chemical_space.get("feature_names", [])
            current_header.extend(self.header)
            self.chemical_space["feature_names"] = current_header

        if len(self.index2rdmol) and self.write_structures:
            self.chemical_space["compounds"] = {}

            for index, rdmol in self.index2rdmol.items():
                if rdmol is not None:
                    self.chemical_space["compounds"][self.index2id[index]] = {self.KEYS.get("smiles", "smiles"): Chem.MolToSmiles(rdmol, True)}

    def _parse_categories(self):
        category2ids = {}
        
        for index, category in self.index2category.items():
            if self.category_field_delimiter is False:
                categories = [category] 
            else:
                categories = [c.strip() for c in category.split(self.category_field_delimiter) if c.strip() != ""]

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

                id2pcp[index] = self._get_pcp_for_rdmol(rdmol)
                i+=1

            empty = [None for x in PROP2LABEL]
            for i, index in enumerate(self.index_order):
                
                if id2pcp.get(index, False):
                    pcps = id2pcp[index]
                else:
                    pcps = empty
                
                self.chemical_space["points"][index][self.KEYS.get("features", "features")].extend(pcps)
                # self.data[i].extend(pcps)                

            current_header = self.chemical_space.get("feature_names", [])
            current_header.extend([PROP2LABEL[prop] for prop in PROPS_ORDER])
            self.chemical_space["feature_names"] = current_header
            # self.original_data = copy.deepcopy(self.data)

    def _get_pcp_for_rdmol(self, rdmol):
        return [round(PROP2FNC[prop](rdmol), 2) for prop in PROPS_ORDER]

    def normalize_data(self, feature_range=(0,1)):
        """Normalizes data to a scale from 0 to 1."""
        print("Data normalization (scale): {}".format(feature_range))

        min_max_scaler = preprocessing.MinMaxScaler(feature_range)
        self.data = min_max_scaler.fit_transform(self.data)
        self.data = [[round(v, 3) for v in row] for row in self.data]

    def _calculate_distance_matrix(self, index_order=None, index2order=None, method=None):
        start = time.time()
        if index_order is None:
            index_order = self.index_order

        if index2order is None:
            index2order = {x: i for i, x in enumerate(self.index_order)}
        
        print("\nCalculating distance matrix: {} compounds".format(len(index_order)))
        
        fps = np.array([self.index2fpobj[index].ToList() for index in index_order])
        self.dist_matrix = metrics.pairwise_distances(fps, metric=self.metric, n_jobs=self.n_jobs)
        end = time.time()
        print("DM SKLEARN:", end - start)

    def _get_edges(self, similarity_threshold=0.7, knn=2, index_order=None, index2order=None, method=None):
        start = time.time()
        self.edges = []
        self.edges_weights = []
        self.index2edges = defaultdict(dict)
        
        if knn is None:
            knn = 0

        if index_order is None:
            index_order = self.index_order

        if index2order is None:
            index2order = {x: i for i, x in enumerate(index_order)}
        
        count = len(index_order)
        print(f"\nCalculating edges [similarity threshold={similarity_threshold}, knn={knn}]: {count} compounds")
        # print("DM INDEXES MATCH", len(index_order) == len(self.dist_matrix))
        self.dist_matrix = 1 - self.dist_matrix
        dims = len(self.dist_matrix)

        if similarity_threshold == 0 and knn == 0:
            print("Get all edges...")
            
            for i, index in enumerate(index_order, 0):
                if i%100 == 0:
                    print("{}/{}".format(i, count))

                if i == dims-1:
                    break

                self.edges.extend([(index2order[index], index2order[idx]) for idx in index_order[i+1:]])
                self.edges_weights.extend([self.dist_matrix[i][j] for j, idx in enumerate(index_order[i+1:], i+1)])

        elif knn == 0:
            for i, index in enumerate(index_order, 0):
                if i%100 == 0:
                    print("{}/{}".format(i, count))

                if i == dims-1:
                    break

                for j, idx in enumerate(index_order[i+1:], i+1):
                    sim = self.dist_matrix[i][j]
                    
                    if sim >= similarity_threshold:
                        self.edges.append((index2order[index], index2order[idx]))
                        self.edges_weights.append(sim)

        else:
            index2edges = defaultdict(list)
            for i, index in enumerate(index_order, 0):

                if i%100 == 0:
                    print("{}/{}".format(i, count))

                if i == dims-1:
                    break

                for j, idx in enumerate(index_order[i+1:], i+1):
                    sim = self.dist_matrix[i][j]
                    
                    if sim >= similarity_threshold:
                        index2edges[index2order[index]].append((index2order[idx], sim))
                        index2edges[index2order[idx]].append((index2order[index], sim))

            for index, edges in index2edges.items():
                edges.sort(key=lambda x: x[1])

                self.edges.extend([(index, x[0]) for x in edges[-knn:]])
                self.edges_weights.extend([x[1] for x in edges[-knn:]])
                    
        end = time.time()
        print("GET EDGES:", end - start)
        print("EDGES: {}".format(len(self.edges)))
        
    def _convert_fps_to_bitvects(self, fps):
        self.index2fpobj = {}
        for i, fp in enumerate(fps):
            self.index2fpobj[i] = self._get_bitvect_for_fp(fp)

    def _get_bitvect_for_fp(self, fp):
        if type(fp) is list and len(fp) == 1:
            fp = fp[0]
        bitvect = cDataStructs.ExplicitBitVect(len(fp))
        on_indexes = [i for i, b in enumerate(fp) if int(b) == 1]
        bitvect.SetBitsFromList(on_indexes)
        return bitvect

    def _pca(self, data, **kwargs):
        pca = decomposition.PCA(n_components=2)
        coords = pca.fit_transform(data)
        return coords

    def _mds(self, data, **kwargs):
        # sklearn implementation
        mds = manifold.MDS(n_components=2, dissimilarity='precomputed')
        coords = mds.fit_transform(data)

        # igraph implementation
        # layout = g.layout_mds(data, 2, arpack_options=igraph.ARPACKOptions(iter=1000))
        # coords = layout.coords
        return coords

    def _fa(self, data, **kwargs):
        fa = decomposition.FactorAnalysis(n_components=2)
        coords = fa.fit_transform(data)
        return coords

    def _isomap(self, data, **kwargs):
        isomap = manifold.Isomap(n_neighbors=200, n_components=2)
        coords = isomap.fit_transform(data)
        return coords

    def _multicore_tsne(self, data, **kwargs):
        tsne = MulticoreTSNE(n_components=2, metric='precomputed', n_jobs=self.n_jobs)
        coords = tsne.fit_transform(data)
        return coords

    def _tsne(self, data, **kwargs):
        tsne = manifold.TSNE(n_components=2, metric='precomputed')
        coords = tsne.fit_transform(data)
        coords = [[float(x[0]), float(x[1])] for x in coords]
        return coords

    def _umap(self, data, **kwargs):
        print("Calculating PCA: 20 components")
        pca = decomposition.PCA(n_components=20)
        data = pca.fit_transform(data)

        # umap = UMAP(n_neighbors=20, min_dist=1, metric="jaccard")
        umap = UMAP(n_neighbors=20, min_dist=1, metric="euclidean")
        coords = umap.fit_transform(data)
        coords = [[float(x[0]), float(x[1])] if not np.isnan(x[0]) else [0, 0] for x in coords]
        return coords

    def _csn(self, data, **kwargs):
        g = igraph.Graph(len(self.index_order))
        print(self.edges, self.edges_weights)
        g.add_edges(self.edges)
        layout = g.layout_fruchterman_reingold(weights=self.edges_weights if kwargs["weights"] else None)            
        coords = layout.coords

        for i, e in enumerate(self.edges):
            self.index2edges[e[0]][e[1]] = round(self.edges_weights[i], 2)

        return coords

    def _csn_scaffolds(self, data, **kwargs):
        g = igraph.Graph(len(self.index_order))
        edges = copy.copy(self.edges)
        edges.extend(self.scaffold_edges)
        g.add_edges(edges)
        
        edges_weights = copy.copy(self.edges_weights)
        edges_weights.extend(self.scaffold_edges_weights)

        layout = g.layout_fruchterman_reingold(weights=edges_weights if kwargs["weights"] else None)
        coords = layout.coords

        for i, e in enumerate(self.edges):
            self.index2edges[e[0]][e[1]] = round(self.edges_weights[i], 2)

        return coords

    def _mst(self, data, **kwargs):
        edges = []

        for i, e in enumerate(self.edges):
            edges.append((e[0], e[1], 1 - self.edges_weights[i]))

        x, y, s, t, _ = tmap.layout_from_edge_list(
            len(data), edges, create_mst=True
        )
        
        coords = list(zip(x, y))
        index2edges = defaultdict(dict)
        
        for i, indexes in enumerate(zip(s, t)):
            self.index2edges[indexes[0]][indexes[1]] = round(self.dist_matrix[indexes[0]][indexes[1]], 2)
        
        return coords

    def _mst_scaffolds(self, data, **kwargs):
        edges = []
        scaffold_index2order = {si: i for i, si in enumerate(self.scaffold_index_order)}
        identity2value = {0: 0.01}
        
        print("SCAFFOLD EDGES")
        for i, e in enumerate(self.edges):
            weight = 1 - self.edges_weights[i]
            edges.append((e[0], e[1], identity2value.get(weight, weight)))

        if not self.only_scaffolds:
            print("SCAFFOLD COMPOUNDS EDGES")
            for i, e in enumerate(self.scaffold_edges):
                edges.append((e[0], e[1], 0))

        x, y, s, t, _ = tmap.layout_from_edge_list(
            len(self.index_order), edges, create_mst=True
        )

        coords = list(zip(x, y))
        index2edges = defaultdict(dict)
        
        for i, ids in enumerate(zip(s, t)):
            sid1 = scaffold_index2order.get(ids[0], False)
            value = None

            if sid1:
                sid2 = scaffold_index2order.get(ids[1], False)

                if sid2:
                    value = round(self.dist_matrix[sid1][sid2], 2)

            self.index2edges[ids[0]][ids[1]] = value

        return coords

    def arrange(self, by="fps", fps=None, method="pca", similarity_threshold=0.7, add_edges=None, weights=False, knn=None, add_scaffolds_category=False, only_scaffolds=False):
        self.dist_matrix = False
        self.edges = []
        self.edges_weights = []
        self.index2edges = defaultdict(dict)
        self.add_scaffolds_category = add_scaffolds_category
        self.only_scaffolds = only_scaffolds

        if knn is not None:
            knn = int(knn)

        similarity_threshold = float(similarity_threshold)
        bitvects = False
        scaffolds_index_order = False

        if type(method) is not list:
            methods = [method]
        else:
            methods = method

            if "," in methods[0]:
                methods = [x.strip() for x in methods[0].split(",")]

        self.add_edges = add_edges
        
        if add_edges in [None, False]:
            self.add_edges = any([METHODS[m].get("edges", False) for m in methods])
        
        for method in methods:
            
            if "scaffolds" in method in  ["csn_scaffolds", "mst_scaffolds"] and scaffolds_index_order is False:
                self.scaffold_index_order, self.scaffold_index2order = self._arrange_by_scaffolds()
                
                try:
                    self._calculate_distance_matrix(
                        index_order=self.scaffold_index_order,
                        index2order=self.scaffold_index2order,
                        method=method
                    )
                except Exception as e:
                    print(e)

                if self.add_edges:
                    self._get_edges(
                        similarity_threshold=similarity_threshold,
                        knn=knn,
                        index_order=self.scaffold_index_order,
                        index2order=self.scaffold_index2order,
                        method=method
                    )

                data = self.data

            elif (by == "dm" or METHODS[method]["dm"]) and self.dist_matrix is False:

                if fps is None:
                    fps = []
                    for index in self.index_order:
                        fps.append(self.index2fpobj[index])
                    
                    bitvects = True
                        
                elif type(fps[0][1]) in [str, int] and not bitvects:
                    self._convert_fps_to_bitvects(fps)
                    bitvects = True
                    
                self._calculate_distance_matrix()

                if self.add_edges:
                    self._get_edges(
                        similarity_threshold=similarity_threshold,
                        knn=knn,
                        method=method
                    )

                # dist_matrix = np.array([np.array(self.dist_matrix[index]) for index in self.index_order])
                data = self.dist_matrix
                by = "dm"

            elif by in ["data", "fps"]:
                if by == "fps":
                    if fps is None and len(self.index2fpobj):
                        fps = []
                        
                        for index in self.index_order:
                            fps.append(self.index2fpobj[index])
                    
                    data = np.array([np.array([int(b) for b in fp]) for fp in fps])
                else:
                    data = np.array([np.array(row) for row in self.data])

            print(f"Calculating {METHODS[method]['label']}...")
            feature_names = [f"{METHODS[method]['dim_label']} {i}" for i in [1, 2]]
            coords = getattr(self, f"_{method}")(data, weights=weights)            

            if method in ["csn", "csn_weighted", "csn_scaffolds", "mst", "mst_scaffolds"] or self.add_edges:
                # if self.dist_matrix is False and self.edges in [False, []]:
                #     self._calculate_distance_matrix(similarity_threshold)

                # if self.edges in [False, []]:
                #     if knn is None:
                #         knn = len(self.index_order)
                #     self._get_edges(similarity_threshold=similarity_threshold, knn=knn)

                for cid, es in self.index2edges.items():
                    if not self.chemical_space["points"][cid].get("links", False):
                        self.chemical_space["points"][cid]["links"] = []

                    for e, weight in es.items():
                        self.chemical_space["points"][cid]["links"].append([e, weight])

            index2coords = {index:coords[i] for i, index in enumerate(self.index_order)}

            for index, values in self.chemical_space["points"].items():
                if index in index2coords:
                    point_features = self.chemical_space["points"][index][self.KEYS.get("features", "features")]
                    features = [round(index2coords[index][0], 5), round(index2coords[index][1], 5)]
                    # features = [index2coords[index][0], index2coords[index][1]]
                    features.extend(point_features)
                    self.chemical_space["points"][index][self.KEYS.get("features", "features")] = features

                else:
                    self.chemical_space["points"].pop(index, None)

            feature_names.extend(self.chemical_space.get("feature_names", []))
            self.chemical_space["feature_names"] = feature_names

    def _arrange_by_scaffolds(self, align_by_scaffold=True):
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
        self.scaffold_edges_weights = []

        for index_1, scaffold in self.index2scaffold.items():
            for index_2 in self.scaffold2index_orders[scaffold]:
                self.scaffold_edges.append((index2order[index_1], index_2))
                self.scaffold_edges_weights.append(1)

        self._add_scaffolds_to_chemical_space()

        return index_order, index2order

    def _add_scaffolds_to_chemical_space(self):
        for index, scaffold in self.index2scaffold.items():
            label = f"Scaffold {index-len(self.data)+1}"

            self.chemical_space["points"][index] = {
                self.KEYS.get("object_ids", "object_ids"): [label],
                self.KEYS.get("label", "label"): label,
                self.KEYS.get("features", "features"): []
            }

            if not self.only_scaffolds:
                self.chemical_space["points"][index]["links"] =  [[x] for x in self.scaffold2indexes[scaffold]]
            else:
                self.chemical_space["points"][index][self.KEYS.get("object_ids", "object_ids")].extend(
                    [str(self.index2id[i]) for i in self.scaffold2index_orders[scaffold]]
                )

            for i, f in enumerate(self.chemical_space.get("feature_names", [])):
                values = [self.chemical_space["points"][x][self.KEYS.get("features", "features")][i] for x in self.scaffold2indexes[scaffold]]
                values = [v for v in values if v is not None]
                value = round(np.mean(values), 2) if len(values) else None
                self.chemical_space["points"][index][self.KEYS.get("features", "features")].append(value)

            self.chemical_space["compounds"][label] = {self.KEYS.get("smiles", "smiles"): scaffold, "color": "red"}
            
            if self.only_scaffolds:
                for i in self.scaffold2index_orders[scaffold]:
                    self.chemical_space["points"].pop(i, None)
        
        if self.only_scaffolds and self.category_field:
            pass
        
        elif self.add_scaffolds_category:
            if not self.chemical_space.get("categories", False):
                self.chemical_space["categories"] = []

            self.chemical_space["categories"].append({
                self.KEYS.get("label", "label"): "scaffolds", 
                "color": "gray", 
                "points": list(self.index2scaffold.keys()), 
                "shape": "circle"
            })


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
            space_json = json.dumps(space_json, ignore_nan=True)
            space_json = self._minify_data(space_json)
        elif dump:
            space_json = json.dumps(space_json, indent=4,  ignore_nan=True)

        if filename:
            output = open(filename, "w")
            output.write(space_json)
        
        return space_json

    def _minify_data(self, data):
        return jsmin.jsmin(data)

def _process_(arguments):
    s = ChemSpace(
        write_structures = False if arguments.dont_write_structures else True,
        fp = arguments.fingerprint,
        category_field = arguments.category_field,
        category_field_delimiter = arguments.category_field_delimiter,
        label_field = arguments.label_field,
        compound_structure_field = arguments.compound_structure_field,
        keep_unparsable_structures = arguments.keep_unparsable_structures,
        fingerprint_field = arguments.fingerprint_field,
        metric = arguments.similarity_metric,
        compressed_data_format=arguments.compressed_data_format,
        round_values=arguments.round_values,
        n_jobs=arguments.n_jobs
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
            add_scaffolds_category=arguments.add_scaffolds_category,
            only_scaffolds=arguments.only_scaffolds
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
    parser.add_argument("-kus", "--keep_unparsable_structures", default=False, help="keep data even if the compound structure is not parsed", action="store_true")
    parser.add_argument("-fp", "--fingerprint", type=str, default="ecfp4", help="fingerprint used for a compound representation (ecfp4, ecfp6, maccs, topological, atom_pairs)")
    # parser.add_argument("-c", "--compounds", type=str, default=None, help="csv(text) compound file with delimited values in a form: id,smiles")
    # parser.add_argument("-cd", "--compounds_delimiter", type=str, default=",", help="delimiter of values in compound file")
    parser.add_argument("-arr", "--arrange_by", default="data", help="arrange data by compound structures (distance matrix) or by input data (data/fps)", type=str)
    parser.add_argument("-asc", "--add_scaffolds_category", default=False, help="add scaffolds as a category (only for arrange by csn_scaffolds)", type=str)
    parser.add_argument("-osc", "--only_scaffolds", default=False, help="add only scaffolds to chemical space, point IDs added as object IDs", action="store_true")
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
    parser.add_argument("-sm", "--similarity_metric", type=str, default="jaccard", help="similarity metric")
    parser.add_argument('-rmc','--remove_columns', nargs='+', default=False, help='columns in data that should not be used')
    parser.add_argument('-cdf','--compressed_data_format', default=False, help='use shorter data keys', action='store_true')
    parser.add_argument("-rv", "--round_values", type=int, default=False, help="the number of decimal places used for rounding")
    parser.add_argument("-njobs", "--n_jobs", type=int, default=1, help="the number of CPUs/threads to use for distance matrix calculation")

    
    args = parser.parse_args()
    _process_(args)


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