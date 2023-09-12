# .bash_profile: Default user startup file for bash shell
#
# Version 1.0, January 15, 2010 (PG)
#
# ---> THIS FILE IS NOT INTENDED TO BE MODIFIED BY THE USER. 
#
# Edit: 
#   ~/.bash.myexport - to define private environment variables, add PATHs, etc.
#   ~/.bash.myalias  - to define private aliases and shell settings
#
###############################################################################
# Set PATH and non-science-related environment variables (shared with sh/ksh)

if [ -e $HOME/.bash.export ]; then
    . $HOME/.bash.export
fi

# Set additional aliases and shell settings (shared with sh/ksh setup)

if [ -e $HOME/.bash.alias ]; then
    . $HOME/.bash.alias
fi

# Set additional science package-related variables/aliases 
# (shared with sh/ksh setup)

if [ -e $HOME/.bash.scienv ]; then
    . $HOME/.bash.scienv
fi

# Clean up environment if duplicate assignments exist

if [ -e ${HOME}/.cleanEnv ]; then
    if [ ! -x ${HOME}/.cleanEnv ]; then
	/bin/chmod u+x ${HOME}/.cleanEnv
    fi
    ${HOME}/.cleanEnv
    if [ -e ${HOME}/.newEnv.sh ]; then
       . ${HOME}/.newEnv.sh
    fi
fi

#################################################################
# Added from https://natelandau.com/my-mac-osx-bash_profile/
export PS1=" ________________________________________________________________\n\w @ \h (\D{%d%b%Y %H:%M:%S}) \n> " 
#export PS1=" ________________________________________________________________\n \w @ \h (\u) \n> "
#export PS1=" _______________________________________________________________________________\n \w @ \h (\u) \n> "
export VISUAL=vim
export EDITOR="$VISUAL"
export CLICOLOR=1
export LSCOLORS=ExFxBxDxCxegedabagacad

alias rm='rm -iv' # ask if overwriting
alias cp='cp -iv' # ask if overwriting
alias mv='mv -iv' # ask if overwriting
alias mkdir='mkdir -pv' #create intermediate dirs if nonexistent
alias ls='ls -Fp' # show directories with trailing /, but not executable *s
alias ll='ls -lAhp' 
alias ..='cd ../'                           # Go back 1 directory level
alias ...='cd ../../'                       # Go back 2 directory levels
alias .3='cd ../../../'                     # Go back 3 directory levels
alias .4='cd ../../../../'                  # Go back 4 directory levels
alias .5='cd ../../../../../'               # Go back 5 directory levels
alias .6='cd ../../../../../../'            # Go back 6 directory levels

alias back="cd -"

defaults write com.apple.LaunchServices LSQuarantine -bool NO

source activate py39
