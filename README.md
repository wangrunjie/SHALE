# SHALE
SHALE: An Efficient Algorithm for Allocation of Guaranteed Display Advertising

Online allocation for Guaranteed Display Advertising

supply:

supply_node     inventory       satisfy_demand

0               40              a

1               40              b

2               40              a,c

3               40              c,b

demand:

demand_node     demand  penalty target_supply

a               40      1.000000        0,2

b               40      1.000000        1,3

c               70      1.000000        2,3

allocation:

demand_node     demand          allocation

a               40              43

b               40              37

c               70              65

remaind:

supply_node     inventory       remaind

0               40              6

1               40              9

2               40              0

3               40              0


Any question mail to wangjie5@xiaomi.com
