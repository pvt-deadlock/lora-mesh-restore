# LoRa Mesh Network Simulation

## License

This software using ns3 (https://gitlab.com/nsnam/ns-3-dev) simulator, which is licensed under the terms of the GNU General Public License v2.0 only (GPL-2.0-only).
See the LICENSE file for more details.

## Table of Contents
* [Building ns-3](#building-ns-3)
* [Running simulation](#running-simulation)

> **NOTE**: Much more substantial information about underlying network simulator about ns-3 can be found at
<https://www.nsnam.org>

## Building ns-3

The code for the framework and the default models provided
by ns-3 is built as a set of libraries. User simulations
are expected to be written as simple programs that make
use of these ns-3 libraries.

To build the set of default libraries and the example
programs included in this package, you need to use the
`ns3` tool. This tool provides a Waf-like API to the
underlying CMake build manager.
Detailed information on how to use `ns3` is included in the
[quick start guide](doc/installation/source/quick-start.rst).

Before building ns-3, you must configure it.
This step allows the configuration of the build options,
such as whether to enable the examples, tests and more.

To configure ns-3 with examples and tests enabled,
run the following command on the ns-3 main directory:

```shell
./ns3 configure --enable-examples --enable-tests
```

Then, build ns-3 by running the following command:

```shell
./ns3 build
```

By default, the build artifacts will be stored in the `build/` directory.

## Running simulation

```shell
./ns3 run scratch/lorawan-advanced-final.cc -- --power=14.0 --sf=10
```

The example contains 2 parameters: power and spreading factor of the interruption. It can also get the random seed with help of "rand" parameter. If there is no seed provided, program will pick random seed based on the current timestamp. 

