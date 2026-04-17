import gdb
import gdb.printing

# Helper function to extract pointer type
def get_data_pointer(val):
    # val['data'] is mfem::Memory<T>
    # val['data']['h_ptr'] is T*
    return val['data']['h_ptr']

# -----------------------------
# mfem::Vector pretty-printer
# -----------------------------
class MfemVectorPrinter:
    def __init__(self, val):
        self.val = val

    def to_string(self):
        size = int(self.val['size'])
        return f"[{size}] mfem::Vector" if size else "empty (mfem::Vector)"

    def children(self):
        size = int(self.val['size'])
        data_ptr = get_data_pointer(self.val)
        elem_type = data_ptr.type.target()
        for i in range(size):
            yield f"[{i}]", (data_ptr + i).dereference()

# -----------------------------
# mfem::Array<T> pretty-printer
# -----------------------------
class MfemArrayPrinter:
    def __init__(self, val):
        self.val = val

    def to_string(self):
        size = int(self.val['size'])
        return f"[{size}] mfem::Array" if size else "empty (mfem::Array)"

    def children(self):
        size = int(self.val['size'])
        data_ptr = get_data_pointer(self.val)
        for i in range(size):
            yield f"[{i}]", (data_ptr + i).dereference()

# -----------------------------
# mfem::DenseMatrix pretty-printer
# -----------------------------
class MfemDenseMatrixPrinter:
    def __init__(self, val):
        self.val = val

    def to_string(self):
        h = int(self.val['height'])
        w = int(self.val['width'])
        return f"({h}, {w}) mfem::DenseMatrix" if w else "empty mfem::DenseMatrix"

    def children(self):
        height = int(self.val['height'])
        width = int(self.val['width'])
        data_ptr = get_data_pointer(self.val)
        for row in range(height):
            for col in range(width):
                idx = row * width + col
                yield f"[{row},{col}]", (data_ptr + idx).dereference()

# -----------------------------
# Register all printers
# -----------------------------
def build_pretty_printer():
    pp = gdb.printing.RegexpCollectionPrettyPrinter("mfem")
    pp.add_printer("mfem::Vector", "^mfem::Vector$", MfemVectorPrinter)
    pp.add_printer("mfem::Array", "^mfem::Array<.*>$", MfemArrayPrinter)
    pp.add_printer("mfem::DenseMatrix", "^mfem::DenseMatrix$", MfemDenseMatrixPrinter)
    return pp

def register_printers(objfile):
    gdb.printing.register_pretty_printer(objfile, build_pretty_printer())
