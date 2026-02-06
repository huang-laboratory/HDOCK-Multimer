import sys


file_models = sys.argv[1]


file_stem, file_suffix = file_models.rsplit(".", 1)
with open(file_models, 'r') as file:
    num_models = 0
    file_current_model = None
    for line in file:
        if line.startswith("MODEL"):
            num_models += 1
            if file_current_model:
                file_current_model.close()
            file_current_model = open(f"{file_stem}_{num_models}.{file_suffix}", 'w')
        elif line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
            if file_current_model:
                print(line.strip(), file=file_current_model)

    if file_current_model:
        file_current_model.close()
