#!/usr/bin/python

from distutils.core import setup

setup(  name='wel',
        version='1.0',
        author='Wilson Freitas',
        author_email='wilson@aboutwilson.net',
        url='http://aboutwilson.net',
        package_dir = {'' : 'Python'},
        packages=['wel'],
        py_modules=['cmz', 'cmzga'],
        data_files=[('data',['data/close.dat']),],
        scripts=['mainsm'])

