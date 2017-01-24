import os
import glob
import shutil

def fix_coltimes(coldir, newdir, resume=None):
    dirs = glob.glob(coldir + '/*/core.*')
    os.makedirs(newdir, mode=0o755, exist_ok=True)
    if resume:
        have_reached_resume = False
    else:
        have_reached_resume = True
    for sim in dirs:
        if not have_reached_resume:
            if sim == resume:
                have_reached_resume = True
            else:
                continue
        hostname = sim
        hostname = hostname[hostname.find('/')+1:]
        hostname = hostname[:hostname.find('/')]
        run_base = os.path.join(newdir, hostname)
        run_num = 0
        while True:
            try:
                run_dir = run_base + '.' + str(run_num)
                os.makedirs(run_dir, mode=0o755)
                break
            except OSError:
                run_num = run_num + 1
        shutil.copytree(os.path.join(sim, 'input'), os.path.join(run_dir, 'input'))
        cpfile = lambda f: shutil.copyfile(os.path.join(sim, f),
                                           os.path.join(run_dir, f))
        cpfile('wlcsim.exe')
        cpfile('run_parameter.pl')
        cpfile('runsim.sh')
        for i,finished_sims in enumerate(glob.glob(sim + '/savedata/*')):
            finished_dir = run_base + '.' + str(run_num + i + 1)
            shutil.copytree(run_dir, finished_dir)
            shutil.copytree(finished_sims, os.path.join(finished_dir, 'data'))
        shutil.copytree(os.path.join(sim, 'data'), os.path.join(run_dir, 'data'))

def fix_old0(old0_dir, newdir):
    dirs = glob.glob(old0_dir + '/tower*/par-run-dir.*/run.*')
    os.makedirs(newdir, mode=0o755, exist_ok=True)
    for sim in dirs:
        if not os.path.isdir(os.path.join(sim, 'input')):
# old.0/tower0/par-run-dir.1/run.25 is a special snowflake
# here is where we would use him if we wanted to, but doesn't seem
# like anything important was saved in this format
            if os.listdir(sim):
                pass
            continue
        hostname = sim
        hostname = hostname[hostname.find('/')+1:]
        hostname = hostname[:hostname.find('/')]
        run_base = os.path.join(newdir, hostname)
        run_num = 0
        while True:
            try:
                run_dir = run_base + '.' + str(run_num)
                os.makedirs(run_dir, mode=0o755)
                break
            except OSError:
                run_num = run_num + 1
        shutil.copytree(os.path.join(sim, 'input'), os.path.join(run_dir, 'input'))
        cpfile = lambda f: shutil.copyfile(os.path.join(sim, f),
                                           os.path.join(run_dir, f))
        try:
            cpfile('wlcsim.exe')
        except FileNotFoundError:
            print('Did not find wlcsim.exe in ' + sim)
        try:
            cpfile('run_parameter.pl')
        except FileNotFoundError:
            print('Did not find run_parameter.pl in ' + sim)
        try:
            cpfile('runsim.sh')
        except FileNotFoundError:
            print('Did not find runsim.sh in ' + sim)
        for i,finished_sims in enumerate(glob.glob(sim + '/savedata/*')):
            finished_dir = run_base + '.' + str(run_num + i + 1)
            shutil.copytree(run_dir, finished_dir)
            shutil.copytree(finished_sims, os.path.join(finished_dir, 'data'))
        shutil.copytree(os.path.join(sim, 'data'), os.path.join(run_dir, 'data'))

