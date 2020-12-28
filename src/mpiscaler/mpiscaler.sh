#!/bin/sh

dirname=$( dirname $0 )
[[ "${dirname:0:1}" != "/" ]] && dirname="$( pwd )/$dirname"
basename=$( basename $0 )

INF=0.00

DATABASE="$dirname/../db/prodb/PRODB"
OUTFILE="$dirname/../db/mpiscaler.I${INF}.`date +'%y%b%d'`.out"
NODEFILE="$dirname/../var/node.lst"
CPUS=4
MPISCALER="$dirname/mpiscaler"

usage="
This is a helper script to launch \`mpiscaler'.
2008(C)Mindaugas Margelevicius,IBT,Vilnius

$basename <Options>

Options:

-d <database>  absolute pathname to database to be read for scaling
       default=$DATABASE
-o <output>    absolute pathname to file to write output to
       default=$OUTFILE
-n <filename>  pathname to file containing list (1/line) of node names
               NOTE: The list should not include this computer
       default=$NODEFILE
-c <number>    number of CPUs to use on nodes
       default=$CPUS
-p <program>   absolute pathname to mpiscaler program
               (must be enabled to be accessed from other nodes)
       default=$MPISCALER
-h             short description

"


while getopts "d:o:n:c:p:h" Option
do
    case $Option in
        d ) DATABASE=${OPTARG} ;;
        o ) OUTFILE=${OPTARG} ;;
        n ) NODEFILE=${OPTARG} ;;
        c ) CPUS=${OPTARG} ;;
        p ) MPISCALER=${OPTARG} ;;
        h ) echo "$usage"; exit 0 ;;
        * ) echo Error: Unrecognized argument.; exit 1 ;;
    esac
done
shift $(( $OPTIND - 1 ))

ERRFILE="$( dirname $OUTFILE )/$( basename $OUTFILE .out ).err"

if [[ -z "$( ls -1 ${DATABASE}* )" ]]; then echo Error: Database does not exist.; exit 1; fi
if [[ -z $NODEFILE || ! -f $NODEFILE ]]; then echo Error: File with the list of nodes does not exist.; exit 1; fi
if [[ -z $MPISCALER || ! -f $MPISCALER ]]; then echo Error: Program $MPISCALER does not exist.; exit 1; fi
if [[ -z $CPUS || $CPUS -le 0 ]]; then echo Error: Wrong number of CPUs specified.; exit 1; fi

NODES=$(< $NODEFILE )


masternodename=$( hostname )
mastercpus=1

##mpiscaler_cmd="$MPISCALER -d $DATABASE -I $INF -o $OUTFILE"
mpiscaler_cmd="$MPISCALER -d $DATABASE -o $OUTFILE"

deployment="-n $mastercpus -host $masternodename $mpiscaler_cmd"



for node in $NODES; do
    [[ "${node:0:1}" == "#" ]] && continue
    nodecpus=$CPUS
    nodecpus=$( echo $node | awk -F ':' '{print $2}')
    nodename=$( echo $node | awk -F ':' '{print $1}')
    [[ -z "$nodecpus" ]] && nodecpus=$CPUS
    deployment="$deployment : -n $nodecpus -host $nodename $mpiscaler_cmd"
done

cmd="mpiexec -l $deployment >$ERRFILE 2>&1 &"

## to run by PBS engine
## cmd="mpiexec -l $deployment"

echo $cmd
eval $cmd

