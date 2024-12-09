import gdb

class GSLMatrixPrinter:
    def __init__(self, val):
        self.val = val

    def children(self):
        rows = int(self.val['size1'])
        cols = int(self.val['size2'])
        tda = int(self.val['tda'])
        data = self.val['data']

        for i in range(rows):
            for j in range(cols):
                index = i * tda + j
                yield f"({i}, {j})", data[index]
    def to_string(self):
        return f"{self.val['size1']} x {self.val['size2']} gsl_matrix_float"

class GSLVectorPrinter:
    def __init__(self, val):
        self.val = val

    def children(self):
        size = int(self.val['size'])
        data = self.val['data']

        for i in range(size):
            yield f"({i})", data[i]
    def to_string(self):
        return f"{self.val['size']} gsl_vector_float"

def register_gsl_printers():
    gdb.pretty_printers.append(lambda val: GSLMatrixPrinter(val) if str(val.type) == "gsl_matrix_float" else None)
    gdb.pretty_printers.append(lambda val: GSLVectorPrinter(val) if (str(val.type) == "gsl_vector_float") or (str(val.type) == "gsl_quat_float") else None)

register_gsl_printers()
