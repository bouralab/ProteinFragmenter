import sys
import numpy as np
import pandas as pd

np_dtypes = {"bigint":int, "varchar":str, "char":str, "int":int,
    "tinyint":int, "double":float, "text":str, "text,":str}

def get_table_properties(f):
    columns = []
    dtypes = {}
    for line in f:
        if line.strip().startswith(")") or line.strip().startswith("KEY") or line.strip().startswith("PRIMARY KEY"):
            break
        fields = line.strip().split()
        column = fields[0][1:-1]
        columns.append(column)
        dtypes[column] = np_dtypes.get(fields[1].split("(",1)[0], object)

    print(columns, dtypes)
    return columns, dtypes

def process_table(info, columns, dtypes, name, chunksize=100):
    # rows = [[dtypes[col](field[1:-1]) if dtypes[col]==str else dtypes[col](field) \
    #     for field, col in zip(row, columns)] for row in \
    #     next(f).split("VALUES", 1)[1].strip()[1:-2].split("),(")]

    rows = []
    for row in info.strip()[1:-2].split("),("):
        fields = row.split(",")
        assert len(fields)==len(columns), "{} {}".format(len(fields),len(columns))
        values = [dtypes[col](field[1:-1]) if dtypes[col]==str else dtypes[col](field) \
            for field, col in zip(fields, columns)]
        rows.append(values)
    print(len(rows), len(columns), len(rows[0]))
    df = pd.DataFrame(rows, columns=columns)
    df.to_hdf("CAPSDB.h5", name, mode="a", append=True, table=True, format="table", complib="bzip2", complevel=9, min_itemsize=1024)
    # print(len(rows))
    # print(rows)
    # assert 0
    # rows = []
    # row = ""
    # clear_first = True
    # while True:
    #     c = f.read(1)
    #
    #     if clear_first and row.endswith("VALUES"):
    #         row = ""
    #         clear_first = False
    #         continue
    #     elif (c == "(" and row[-2:] == "),") or c==";":
    #         fields = row[1:-2].split(",")
    #         assert len(fields)==len(columns)
    #         fields = [dtypes[col](field[1:-1]) if dtypes[col]==str else dtypes[col](field) \
    #             for field, col in zip(fields, columns)]
    #         if len(rows) < chunksize:
    #             rows.append(fields)
    #         else:
    #             df = pd.DataFrame(rows, columns=columns)
    #             print(df)
    #             df.to_hdf("CAPSDB.h5", name, mode="a", table=True, format="table", complib="bzip2", complevel=9, min_itemsize=1024)
    #             del df
    #             del rows
    #             rows = []
    #         row = "("
    #         if c==";":
    #             break
    #         continue
    #     else:
    #         row += c
    # print("OUT")
    # if rows:
    #     print("SAVE LAST")
    #     df = pd.DataFrame(rows, columns=columns)
    #     df.to_hdf("CAPSDB.h5", name, mode="a", table=True, format="table", complib="bzip2", complevel=9, min_itemsize=1024)
    #     del df
    #     del rows
    #     rows = []

def process_caps_db(capsdb_file):
    allcolumns = {}
    with open(capsdb_file) as f:
        columns, dtypes, name = None, None, None
        for line in f:
            if line.startswith("CREATE TABLE"):
                name = line.split()[2][1:-1]
                columns, dtypes = get_table_properties(f)
                allcolumns[name] = (columns, dtypes)
            if line.startswith("INSERT INTO"):
                header, info = line.split("VALUES", 1)
                name = header.strip().split()[2][1:-1]
                columns, dtypes = allcolumns[name]
                process_table(info, columns, dtypes, name)

if __name__ == "__main__":
    process_caps_db(sys.argv[1])
