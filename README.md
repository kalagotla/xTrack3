----------------------------------------------------------------------
#7                                                                  7#
#7                            xTrack3                               7#
#7           Experimental particle Tracker based on VISUAL3         7#
#7        By : Dilip Kalagotla ~ kal @ dilip.kalagotla@gmail.com    7#
#7                         Date: 04-15/2021                         7#
----------------------------------------------------------------------

File Structure:
xTrack3
├── build
├── lib
│   ├── inc
│   │   ├── getblock.inc
│   │   ├── Visual3.inc
│   │   └── XFtn.inc
│   ├── libVisual3.a
│   └── libVisual3OLD.a
├── Makefile
├── output
│   └── umData
│       └── 281
├── README.md
├── src
│   └── streamLines
│       └── StreamLines.f
└── test
    └── umPIV
        ├── data
        │   ├── flowData.txt
        │   ├── gridData.txt
        │   └── spec.col
        ├── flowData.f90
        ├── gridData.f90
        ├── init.f
        └── Makefile

Makefile:
	xTrack3 --> make; will create libVisual3.a file
	umPIV   --> make; will create init executable to run UM-PIV case

StreamLines.f --> used to modify particle physics and generator mechanism
init.f        --> used to modify the program behavior


