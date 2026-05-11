#!/bin/sh

# A filter for use like this:
#
#
#  ./describe_demo.sh | ./filter.sh
#
# to make output like this:
#
# (venv) CGL-M4-2222:symcoder lester$ ./describe_demo.sh | ./filter.sh
# [0:3]    ORBIT      mag   .     .   u=(1)  v=(.)  shared=(.)  len=3  .                |  mag(a)
# [3:6]    ORBIT      dot   .     .   u=(2)  v=(.)  shared=(.)  len=3  .                |  dot(a,b)
# [6:7]    ORBIT      eps3  .     .   u=(3)  v=(.)  shared=(.)  len=1  sign_compressed  |  eps3(a,b,c)
# [7:7]    NULL_COMP  dot   dot   SS  u=(2)  v=(2)  shared=(1)  len=0  .                |  dot(a,b)×dot(a,c)
# [7:7]    NULL_SELF  dot   dot   SS  u=(2)  v=(2)  shared=(2)  len=0  .                |  dot(a,b)×dot(a,b)
# [7:7]    NULL_COMP  dot   eps3  SA  u=(2)  v=(3)  shared=(2)  len=0  .                |  dot(a,b)×eps3(a,b,c)
# [7:13]   ASSOC      dot   mag   SS  u=(2)  v=(1)  shared=(0)  len=6  .                |  dot(a,b)×mag(c)
# [13:13]  NULL_COMP  dot   mag   SS  u=(2)  v=(1)  shared=(1)  len=0  .                |  dot(a,b)×mag(a)
# [13:13]  NULL_SELF  eps3  eps3  AA  u=(3)  v=(3)  shared=(3)  len=0  .                |  eps3(a,b,c)×eps3(a,b,c)
# [13:13]  NULL_COMP  eps3  mag   AS  u=(3)  v=(1)  shared=(1)  len=0  .                |  eps3(a,b,c)×mag(a)
# [13:13]  NULL_COMP  mag   mag   SS  u=(1)  v=(1)  shared=(0)  len=0  .                |  mag(a)×mag(b)
# [13:13]  NULL_SELF  mag   mag   SS  u=(1)  v=(1)  shared=(1)  len=0  .                |  mag(a)×mag(a)


grep -A 100000 START | grep -B 100000 STOP | grep -v === | column -t


