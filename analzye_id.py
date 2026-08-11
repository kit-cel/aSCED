'''
Script to transform 
'''
def binary_row_to_hex(binary_row):
    """
    Convert a binary row (list of 0/1) to a hexadecimal string.
    """
    bit_string = "".join(map(str, binary_row))
    return hex(int(bit_string, 2))[2:].upper()


def hex_to_binary_row(hex_identifier, row_length):
    """
    Invert a hexadecimal string back to a binary row of given length.
    """
    # Convert hex to integer, then to binary without '0b' prefix
    bit_string = bin(int(hex_identifier, 16))[2:]

    # Restore leading zeros if necessary
    bit_string = bit_string.zfill(row_length)

    return list(map(int, bit_string))


# -----------------------
# Test / Demonstration for BCH 63,36
# -----------------------
if __name__ == "__main__":
    # Example binary row (may include leading zeros)
    original_row = [0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0]

    ids = [
        "3D5F09BC4506449B",
        "79EBB62DE5975465",
        "134738CE59939DD4",
        "4F72364706BDB4E1",
        "4F72364706BDB4E1",
        "435D604050EC27C5",
        "27CAF3F27ED00830",
        "54F9129F7DA9A852",
        "66EC41F0CAA9DB96",
        "60BBD2CA3082C079",
        "4F72364706BDB4E1",
        "79EBB62DE5975465",
        "435D604050EC27C5",
        "48FC93BCEDBFC4BA",
        "3D5F09BC4506449B",
    ]

    for gid in ids:
        recovered = hex_to_binary_row(gid, 63)
        print(f"Good ID: {gid} -> Recovered binary row: {recovered}")
        print("Row weight:", sum(recovered))

    # Forward conversion
    # hex_id = binary_row_to_hex(original_row)

    # # Inverse conversion
    # recovered_row = hex_to_binary_row(hex_id, len(original_row))

    # print("Original binary row: ", original_row)
    # print("Hex identifier:      ", hex_id)
    # print("Recovered binary row:", recovered_row)
    # print("Row weight:", sum(recovered_row))

    # # Verification
    # assert original_row == recovered_row
    # print("Test passed: inversion is correct.")


