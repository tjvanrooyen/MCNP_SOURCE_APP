@echo off
setlocal EnableExtensions
title Commit to GitHub
REM ===========================================================================
REM  Commit_to_GitHub.cmd
REM
REM  Commits ALL changes in this project folder and pushes them to GitHub.
REM  Double-click it, or run it from a terminal.  It:
REM    1. shows what has changed;
REM    2. asks for a commit message  [Enter = a dated default message];
REM    3. commits everything  [files listed in .gitignore are skipped];
REM    4. first pulls in any changes made on GitHub itself, then pushes.
REM
REM  On first use it asks for the name and e-mail to record on commits and
REM  stores them for this repository only.  The first push may open a GitHub
REM  sign-in window.
REM ===========================================================================

cd /d "%~dp0"

where git >nul 2>&1
if errorlevel 1 goto :no_git
git rev-parse --is-inside-work-tree >nul 2>&1
if errorlevel 1 goto :no_repo

for /f "delims=" %%B in ('git rev-parse --abbrev-ref HEAD') do set "BRANCH=%%B"

REM --- Author identity: asked once, then stored for this repository only ----
set "GIT_NAME="
set "GIT_EMAIL="
for /f "delims=" %%A in ('git config user.name 2^>nul') do set "GIT_NAME=%%A"
for /f "delims=" %%A in ('git config user.email 2^>nul') do set "GIT_EMAIL=%%A"
if not defined GIT_NAME call :ask_name
if not defined GIT_EMAIL call :ask_email
if not defined GIT_NAME goto :no_identity
if not defined GIT_EMAIL goto :no_identity

echo.
echo ===========================================================================
echo  Repository : %CD%
echo  Branch     : %BRANCH%  --^>  GitHub
echo  Author     : %GIT_NAME% ^<%GIT_EMAIL%^>
echo ===========================================================================
echo.

REM --- Commit -----------------------------------------------------------------
set "CHANGES="
for /f "delims=" %%S in ('git status --porcelain') do set "CHANGES=1"
if not defined CHANGES (
    echo No local changes to commit.
    goto :sync
)

echo Changes to be committed:
git status --short
echo.
set "MSG="
set /p "MSG=Commit message [press Enter for a dated default]: "
if not defined MSG for /f "delims=" %%D in ('powershell -NoProfile -Command "Get-Date -Format 'yyyy-MM-dd HH:mm'"') do set "MSG=Update %%D"

REM Write the message to a file so quotes, ampersands etc. survive intact.
set "MSGFILE=%TEMP%\Commit_to_GitHub_msg.txt"
setlocal EnableDelayedExpansion
> "!MSGFILE!" echo(!MSG!
endlocal

git add -A
if errorlevel 1 goto :fail
git commit -q -F "%MSGFILE%"
if errorlevel 1 goto :fail
del "%MSGFILE%" >nul 2>&1
echo.
echo Committed:
git log --oneline -1

REM --- Bring in changes made on GitHub, then push ----------------------------
:sync
echo.
git ls-remote --exit-code --heads origin %BRANCH% >nul 2>&1
if errorlevel 1 goto :push
echo Fetching any changes made on GitHub ...
git -c advice.mergeConflict=false pull --rebase origin %BRANCH%
if errorlevel 1 goto :pull_failed

:push
echo.
echo Pushing to GitHub ...
git push -u origin %BRANCH%
if errorlevel 1 goto :push_failed

echo.
echo Done.  Latest commits on GitHub:
git log --oneline -3
echo.
pause
exit /b 0

REM --- Helpers and error exits ------------------------------------------------
:ask_name
set /p "GIT_NAME=Your name for commits: "
if defined GIT_NAME git config user.name "%GIT_NAME%"
exit /b

:ask_email
set /p "GIT_EMAIL=Your e-mail for commits: "
if defined GIT_EMAIL git config user.email "%GIT_EMAIL%"
exit /b

:no_git
echo ERROR: git was not found.  Install Git for Windows and try again.
goto :fail

:no_repo
echo ERROR: this folder is not a git repository.
goto :fail

:no_identity
echo ERROR: a name and an e-mail are needed to record who made the commit.
goto :fail

:pull_failed
git rebase --abort >nul 2>&1
echo.
echo ERROR: could not bring in the changes made on GitHub.  Either GitHub could
echo        not be reached, or the same lines were edited both here and on
echo        GitHub.  Nothing was pushed; your commit is safe on this PC.
goto :fail

:push_failed
echo.
echo ERROR: the push to GitHub failed - see the message above.
echo        Your commit is safe on this PC; run this script again to retry.
goto :fail

:fail
echo.
pause
exit /b 1
