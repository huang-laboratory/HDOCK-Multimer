from argparse import ArgumentParser
from itertools import combinations_with_replacement
from pathlib import Path, PurePath
import os
import shutil
import subprocess
from subprocess import run, PIPE
from concurrent.futures import ProcessPoolExecutor
import gc
from tempfile import TemporaryDirectory
from datetime import datetime
import time
import multi_utils
from picker import Picker
from chainlist import ChainList


SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))
scorecom = os.path.join(TOOL_PATH, "scorecom.sh")
hdock = os.path.join(TOOL_PATH, "hdock")
createpl = os.path.join(TOOL_PATH, "createpl")

current_dir = os.getcwd()


def process(cmd: str):
    run(cmd, shell=True, check=True)

def is_relative_to(path, base):
    try:
        Path(path).relative_to(base)
        return True 
    except ValueError:
        return False

def get_num_cpus(max_workers: int):
    cpu_cap = int(os.getenv("NUM_CPUS", max_workers))
    cpu_count = os.cpu_count()
    return max(1, min(cpu_count, cpu_cap))


class MultiDock:
    """
    Multibody docking protocol,
    all properties are initialized by command line arguments
    """
    def __init__(self, fout) -> None:
        args = self.arg_parser()
        print()
        print(args)

        self.fout = fout
        self.rounds = args.rounds
        self.intermediates = args.intermediates
        self.no_itscore = args.no_itscore
        self.no_append = args.no_append
        self.proc = args.proc
        self.distri = args.distri
        self.base: str = args.base
        self.num_cpus = get_num_cpus(self.proc)
        self.picker = Picker.get_picker(
            self.fout, args.norm, args.pick, args.model, not args.no_append
        )
        self.__chains = ChainList(args.chains, args.repeat)
        self.chain_num = len(args.chains)
        if self.rounds is None:
            self.rounds = len(self.__chains)

        self.cwd = Path.cwd()
        # proc = run("hostname", stdout=PIPE, check=True)
        # self.hostname = proc.stdout.decode()[:-1]

    def arg_parser(self):
        """ return arguments namespace """
        arg_parser = ArgumentParser()
        arg_parser.add_argument(
            "chains", nargs="+", help="pdb files of each chain", type=multi_utils.pdb_file_checker
        )
        arg_parser.add_argument(
            "-i", "--intermediates", help="retrieve all intermediates", action="store_true"
        )
        arg_parser.add_argument(
            "-n", "--norm", help="normalize scores by interface atoms", action="store_true"
        )
        arg_parser.add_argument(
            "--no-append", help="cannot append in assembling", dest="no_append", action="store_true"
        )
        arg_parser.add_argument(
            "--no-itscore", help="not use IT-score in HDOCK", dest="no_itscore", action="store_true"
        )
        arg_parser.add_argument(
            "-p", "--pick", type=int, choices=[0, 1, 2, 3], default=2, help="pick model id (default:%(default)s)"
        )
        arg_parser.add_argument(
            "-s", type=int, nargs="*", dest="repeat",
            help="stoichiometry, values are corresponding to chains respectively. If not set, each chain will be 1"
        )
        arg_parser.add_argument(
            "-c", "--complex-num", help="number of complex", type=int, dest="complex_num"
        )
        arg_parser.add_argument(
            "-m", "--model", help="number of models", type=int, choices=range(1, 11), default=1, metavar="1-10"
        )
        arg_parser.add_argument(
            "-r", "--rounds", type=int,
            help="round number of pairwise docking, program may stop earlier if all chains are docked"
        )
        arg_parser.add_argument(
            "-j", "--proc", type=int, default=1, help="maximal process number"
        )
        arg_parser.add_argument(
            "-d", "--distri", dest="distri", action="store_true", default=False,
            help="distribute docking tasks to multiple nodes"
        )
        arg_parser.add_argument(
            "-b", "--base", dest="base", type=str, default=None, help="hdock pairwise output folder to reuse"
        )
        args = arg_parser.parse_args()

        if args.repeat is None:
            if len(args.chains) < 2:
                arg_parser.error("[ERROR] Number of chains must be greater than 1.")
        elif len(args.repeat) != len(args.chains):
            arg_parser.error("[ERROR] Stoichiometry is not matched with chains.")

        if args.complex_num is None:
            args.complex_num = len(args.chains)
        elif args.complex_num > len(args.chains):
            arg_parser.error("[ERROR] Complex number is larger than number of chains.")
        return args

    def __dock_cmds(self):
        """ pairwise docking command generator """
        itscore = " -itscore false" if self.no_itscore else ""

        for rec, lig in combinations_with_replacement(self.__chains.init_chains, 2):
            if rec == lig and self.__chains.chain_repeat(rec) == 1:
                continue
            complex_name = multi_utils.get_outfile_name(rec, lig)
            if Path(complex_name).exists():
                continue

            if (
                self.base is not None
                and Path(self.cwd / self.base / complex_name).exists()
            ):
                shutil.copy(self.cwd / self.base / complex_name, complex_name)
                yield (
                    f"{createpl} {complex_name} {multi_utils.get_model_file(rec, lig, 1)} -nmax 1"
                )
            else:
                # node = "" if self.distri else f"-w {self.hostname}"
                yield (
                    f"{hdock} {rec} {lig} -out {complex_name}{itscore} > {complex_name}_hdock && "
                    f"{createpl} {complex_name} {multi_utils.get_model_file(rec, lig, 1)} -nmax 1"
                )

    def __pairwise(self):
        """ run pairwise docking and return a list of new files for next docking round """
        # sort chains by total atoms in pdb, facilitating better docking performance
        print(self.__chains, file=self.fout)

        stime = datetime.now()

        dock_cmds = list(self.__dock_cmds())
        proc = int(min(len(dock_cmds), self.num_cpus))
        # if proc > 0:
        #     with multiprocessing.Pool(proc) as pool:
        #         pool.map(process, dock_cmds)
        if proc > 0:
            with ProcessPoolExecutor(max_workers=proc) as executor:
                futures = {executor.submit(process, cmd): cmd for cmd in dock_cmds}
            for future in futures:
                try:
                    future.result()
                except Exception as e:
                    print(f"[ERROR] Error in executing command: {futures[future]}", flush=True)
                    print(f"[ERROR] Exception: {e}", flush=True)
            gc.collect()

        etime = datetime.now()
        print("Dock in this round consumed: ", etime - stime, file=self.fout)

        # pick chains which can form a complex
        return self.picker.pick(self.__chains)

    def start(self):
        """ start multibody docking until all chains have been bound to one complex """
        # print(f"[INFO] {self.hostname} node", file=self.fout)
        n_round = 0
        if len(self.__chains) > 1 and self.rounds > 0:
            start_time = time.time()
            with TemporaryDirectory() as tempdir:
                # copy chains to temporary dir and remain the relative path
                for i in self.__chains.init_chains:
                    os.makedirs(tempdir / PurePath(i).parent, exist_ok=True)
                    shutil.copy(self.cwd / i, tempdir / PurePath(i).parent)

                os.chdir(tempdir)
                Path(multi_utils.DEFAULT_MODEL_FOLDER).mkdir(parents=True, exist_ok=True)
                self.__chains.duplicate()

                local_model_path = self.cwd / multi_utils.DEFAULT_MODEL_FOLDER
                local_model_path.mkdir(parents=True, exist_ok=True)
                for _ in range(self.rounds):
                    n_round += 1
                    print("\n************************************************************", file=self.fout)
                    print(f"current round: {n_round}", file=self.fout)
                    print("chains to be docked:", file=self.fout)
                    self.__chains = self.__pairwise()
                    self.__chains.duplicate()

                    # copy all new chains to local after each pairwise docking round
                    if self.intermediates:
                        print("retrieve all intermediates", file=self.fout)
                        for chain in self.__chains:
                            #if Path(chain).is_relative_to(multi_utils.DEFAULT_MODEL_FOLDER):
                            if is_relative_to(chain, multi_utils.DEFAULT_MODEL_FOLDER):
                                if len(self.__chains) == 1:
                                    shutil.copy(chain, local_model_path / "multi.pdb")
                                else:
                                    shutil.copy(chain, local_model_path)

                    if len(self.__chains) == 1:
                        if not self.intermediates:
                            shutil.copy(self.__chains[0], local_model_path / "multi.pdb")
                        break

                # copy all pairwise docking output files to local
                if self.intermediates:
                    local_outfile_path = self.cwd / multi_utils.DEFAULT_OUTFILE_PATH
                    local_outfile_path.mkdir(parents=True, exist_ok=True)
                    for outfile in Path(".").glob("*.out"):
                        shutil.copy(outfile, local_outfile_path)

            print(self.__chains)
            final_model_path = local_model_path / "multi.pdb"
            result = run(f"bash {scorecom} {final_model_path}", shell=True, check=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            it_score = float(result.stdout.strip())
            end_time = time.time()
            print(f"\nAsymmetric docking cost {end_time - start_time}", flush=True)
            print(f"\nMulti-body docking cost {end_time - start_time}",
                  file=self.fout, flush=True)
            print(f"\nIT-score: {it_score}", file=self.fout, flush=True)


if __name__ == "__main__":
    with open(f"{current_dir}/multi.out", 'w+') as fout:
        MultiDock(fout).start()
