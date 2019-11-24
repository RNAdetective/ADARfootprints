#!/usr/bin/env bash
create_dir() {
if [ ! -d "$new_dir" ];
then
  mkdir "$new_dir"
fi
}
install_tool() {
sudo apt-get --yes install "$intool"
}
checktools() {
if ! [ -x "$(command -v muscle)" ]; then
  echo 'Error: muscle is not installed.' >&2
  sudo apt-get install muscle
fi
if ! [ -x "$(command -v transeq)" ]; then
  echo 'Error: emboss is not installed.' >&2
  sudo apt-get install emboss
fi
if ! [ -x "$(command -v csvtool)" ]; then
  echo 'Error: emboss is not installed.' >&2
  sudo apt-get install csvtool
fi
}
temp_file() {
if [ -s "$home_dir"/temp.csv ];
then
  rm "$file_in"
  mv "$home_dir"/temp.csv "$file_in"
fi
}
echo "Would you like to setup ADARfootprints on your ubuntu 18? yes or no"
read answer
echo "Please enter your home directory example /home/user"
read home_dir
if [ "$answer" == "yes" ];
then
  sudo apt-get --yes update
  sudo apt-get --yes upgrade
  checktools
  file_in="$home_dir"/ADARfootprints/ADARfootprints.desktop
  sed 's/home_dir/'$home_dir'/g' "$file_in" >> "home_dir"/temp.desktop
  temp_file
  cp "$file_in" "$home_dir"/Desktop
  chmod +x "$home_dir"/Desktop/ADARfootprints.desktop
  gsettings set org.gnome.desktop.background picture-uri "file://"$home_dir"/ADARfootprints/ADARfootprints.jpg"
fi
  
