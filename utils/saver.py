import pickle
import json
import dill

def save_via_pickle(something, name):
    name = "../save_data/" + name + ".pkl"
    with open(name, "wb") as f:
        pickle.dump(something, f)


def load_via_pickle(name):
    name = "../save_data/" + name + ".pkl"
    with open(name, "rb") as f:
        return pickle.load(f)

def save_via_json(something, name):
    name = "../save_data/" + name + ".json"
    with open(name, "wb") as f:
        json.dump(something, f)

def load_via_json(name):
    name = "../save_data/" + name + ".json"
    with open(name, "rb") as f:
        return json.load(f)

def save_via_dill(something, name):
    name = "../save_data/" + name + ".dill"
    with open(name, "wb") as f:
        dill.dump(something, f)

def load_via_dill(name):
    name = "../save_data/" + name + ".dill"
    with open(name, "rb") as f:
        return dill.load(f)

def save(something, name):
    save_via_dill(something, name)

def load(name):
    load_via_dill(name)