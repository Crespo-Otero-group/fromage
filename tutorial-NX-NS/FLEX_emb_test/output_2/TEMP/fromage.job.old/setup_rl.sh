total=$(($(head -1 mol.init.xyz) + $(head -1 shell.xyz)))
echo $total > rl/rl.temp
echo "shell.xyz" >> rl/rl.temp
echo "    XXX__POS__XXX" >> rl/rl.temp
tail -n +3 shell_fixed.xyz >> rl/rl.temp
