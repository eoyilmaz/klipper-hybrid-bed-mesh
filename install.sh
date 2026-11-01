#!/bin/bash
# Links hybrid_bed_mesh.py into Klipper and optionally adds a Moonraker
# update_manager section, "./install.sh -u" removes both.
maintainer=eoyilmaz
repo=klipper-hybrid-bed-mesh
module_name=hybrid_bed_mesh.py

script_path=$(realpath "$0")
repo_path=$(dirname "$script_path")

klipper_path=${KLIPPER_PATH:-~/klipper}
extras_path="$klipper_path/klippy/extras"
moonraker_conf=${MOONRAKER_CONF:-~/printer_data/config/moonraker.conf}

# --------------------------------------------------------------
# Check if running as root
if [ "$(id -u)" = "0" ]; then
    echo "Script must run from non-root !!!"
    exit 1
fi

if [ ! -d "$extras_path" ]; then
    echo "Klipper not found at $klipper_path, set KLIPPER_PATH to your Klipper folder."
    exit 1
fi

# --------------------------------------------------------------
# Uninstall
if [ "$1" = "-u" ]; then
    if [ -L "$extras_path/$module_name" ]; then
        rm "$extras_path/$module_name"
        echo "Removed the $module_name link from $extras_path"
    fi
    if [ -f "$moonraker_conf" ] && grep -q "^\[update_manager $repo\]$" "$moonraker_conf"; then
        sudo service moonraker stop
        # the section runs until the next blank line
        sed -i "/^\[update_manager $repo\]$/,/^\s*$/d" "$moonraker_conf"
        sudo service moonraker start
        echo "Removed [update_manager $repo] from $moonraker_conf"
    fi
    echo "Remove the [hybrid_bed_mesh] section from your printer.cfg and restart Klipper."
    exit 0
fi

# --------------------------------------------------------------
# Link the Klipper module, always force the symlink
ln -sf "$repo_path/extras/$module_name" "$extras_path/"
echo "Linked $module_name to $extras_path"

# --------------------------------------------------------------
# Add the update block to moonraker.conf
if [ -f "$moonraker_conf" ]; then
    if ! grep -q "^\[update_manager $repo\]$" "$moonraker_conf"; then
        read -r -p " Do you want to install the updater? (y/n): " answer
        if [ "$answer" != "${answer#[Yy]}" ]; then
            sudo service moonraker stop
            {
                echo ""
                echo "[update_manager $repo]"
                echo "type: git_repo"
                echo "path: $repo_path"
                echo "origin: https://github.com/$maintainer/$repo.git"
                echo "primary_branch: main"
                echo "managed_services: klipper"
            } >> "$moonraker_conf"
            echo "Added [update_manager $repo] to $moonraker_conf"
            sudo service moonraker start
        else
            echo "Installing the updater skipped"
        fi
    else
        echo "[update_manager $repo] already exists in $moonraker_conf"
    fi
fi

echo "Add a [hybrid_bed_mesh] section to your printer.cfg (see config/example.cfg) and restart Klipper."
