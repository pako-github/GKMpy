#!/bin/python3

class EdgeDoesNotExists(Exception):
    """EdgeDoesNotExists"""
    def __init__(self, msg):
        super(EdgeDoesNotExists, self).__init__(msg)

class NoCompatibleConnectionFound(Exception):
    """EdgeDoesNotExists"""
    def __init__(self, msg):
        super(NoCompatibleConnectionFound, self).__init__(msg)

