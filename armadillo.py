#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    # Find armadillo
    if conf.options.armadillo_dir:
        if not conf.options.armadillo_incdir:
            conf.options.armadillo_incdir=conf.options.armadillo_dir + "/include"
        if not conf.options.armadillo_libdir:
            conf.options.armadillo_libdir=conf.options.armadillo_dir + "/lib"

    armadillo_incdir=[]
    if conf.options.armadillo_incdir:
        armadillo_incdir=[conf.options.armadillo_incdir]
        
    armadillo_libdir=[]
    if conf.options.armadillo_libdir:
        armadillo_libdir=[conf.options.armadillo_libdir]

    armadillo_libs=['armadillo']
    if conf.options.armadillo_libs:
        armadillo_libs=conf.options.armadillo_libs.split()

    conf.check_cxx(msg="Checking for Armadillo",
                          header_name='armadillo',
                          includes=armadillo_incdir,
                          uselib_store='armadillo',
                          libpath=armadillo_libdir,
                          rpath=armadillo_libdir,
                          lib=armadillo_libs)
        
def options(opt):
    armadillo=opt.add_option_group('Armadillo Options')
    armadillo.add_option('--armadillo-dir',
                   help='Base directory where armadillo is installed')
    armadillo.add_option('--armadillo-incdir',
                   help='Directory where armadillo include files are installed')
    armadillo.add_option('--armadillo-libdir',
                   help='Directory where armadillo library files are installed')
    armadillo.add_option('--armadillo-libs',
                   help='Names of the armadillo libraries without prefix or suffix\n'
                   '(e.g. "armadillo")')
