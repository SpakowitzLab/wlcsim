import os

def path_parts(path):
    folders = []
    while True:
        path, folder = os.path.split(path)
        if folder != "":
            folders.append(folder)
        else:
            if path != "":
                folders.append(path)
            break
    folders.reverse()

def get_unique_folder(basename):
    num = 0
    while True:
        try:
            folder = basename + '.' + str(num)
            os.makedirs(folder, mode=0o755)
            break
        except OSError:
            num = num + 1
    return folder
