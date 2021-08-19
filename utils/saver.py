import pickle

def save(something, name):
    name = "../save_data/" + name + ".pkl"
    with open(name, "wb") as f:
        pickle.dump(something, f)


def load(name):
    name = "../save_data/" + name + ".pkl"
    with open(name, "rb") as f:
        return pickle.load(f)
