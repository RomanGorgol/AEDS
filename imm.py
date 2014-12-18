import os
import re
from fabric.api import abort, cd, env, get, hide, hosts, local, prompt, \
    put, require, roles, run, runs_once, settings, show, sudo, warn, puts

from fabric.colors import green, blue, cyan, magenta, red, white, yellow


__author__ = 'pahaz'


def local_env():
    env.user = 'vagrant'
    env.password = 'vagrant'
    env.host_string = '127.0.0.1:2200'


def supercomputer_env():
    from imm_config import PASSWORD, USER, SERVER
    env.user = USER
    env.password = PASSWORD
    env.host_string = SERVER
    # um64.imm.uran.ru
    # -> umu30


def _put(name):
    if os.path.isdir(name):
        run('mkdir -p ~/' + name)
    put(name, '~/', use_sudo=False)


def _get(name):
    get('~/' + name, '.')


def _mpi(name):
    """MPI with OpenMP compile

    Example:

       _mpi('mpi.c')

    """
    if os.path.extsep not in name:
        raise Exception('filename does not contain extsep')
    base_name, name_ext = name.rsplit(os.path.extsep, 1)

    _put(name)
    run('mpicc {name} -openmp -o {base_name} && ./{base_name}'
        .format(**locals()))

    run('mqrun -np 2 -stdout out.txt -stderr err.txt ./{base_name}')
    run('tail ./err.txt')
    run('tail ./out.txt')
    _get('out.txt')
    _get('err.txt')


def _cuda(name):
    """Cuda compile

    Example:

        _cuda('cuda.cu')

    """
    if os.path.extsep not in name:
        raise Exception('filename does not contain extsep')
    base_name, name_ext = name.rsplit(os.path.extsep, 1)

    _put(name)
    # run('nvcc {name} -o {base_name} && ./{base_name}'.format(**locals())
    run('nvcc {name} -o {base_name} && srun --gres=gpu:2 ./{base_name}'
        .format(**locals()))


def main_example():
    run('rm -r kmeans || echo "No!"')
    _put('kmeans')
    C = 10
    POINTS = 1000
    run('g++ ./kmeans/data-gen.cpp -o ./kmeans/data-gen && ./kmeans/data-gen 1 {0} {1} ./kmeans/data-gen.txt'.format(C, POINTS))
    run('g++ ./kmeans/kmeans-parallel.cpp -o ./kmeans/kmeans.elf -fopenmp && ./kmeans/kmeans.elf {0} ./kmeans/data-gen.txt ./kmeans/kmeans.txt'.format(C))
    _get('kmeans')


def main():
    _put('main.c')
    _put('solver.c')
    _put('config.h')


if __name__ == "__main__":
    supercomputer_env()
    main()
