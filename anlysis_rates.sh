g++ -std=c++17 /home/espedica/2025tb/WiP_v_1_1_x/2_stability_plot.cpp -o stab $(root-config --cflags --libs) -lstdc++fs
g++ -std=c++17 /home/espedica/2025tb/WiP_v_1_1_x/control_plots.cpp -o cont $(root-config --cflags --libs) -lstdc++fs
#g++ -std=c++17 /home/espedica/2025tb/WiP_v_1_1_x/nocut5mrad_control_plots.cpp -o nocut5mrad_cont $(root-config --cflags --libs) -lstdc++fs

for i in $@; do
 mkdir -p "plots/run$i"
 ./stab "txt/run$i/" "$i"
 ./cont "txt/run$i/" "$i" 0
 #./nocut5mrad_cont "txt/run$i/" "$i" 0
 echo ">>> Run $i completato"
 echo
done
