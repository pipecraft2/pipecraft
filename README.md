# PipeCraft2 <img src='src/assets/PipeCraft2_logo.png' align="right" height="180" />

## Release 1.0.0

Download: https://github.com/pipecraft2/pipecraft/releases/tag/v1.0.0

---

**USER GUIDE**: https://pipecraft2-manual.readthedocs.io/en/latest/


---

## For developers

Dockerhub: https://hub.docker.com/u/pipecraft

### Prerequisites:

NodeJS https://nodejs.org/en/download/ (make sure you lets node install build tools on windows and build-essential on ubuntu)
Yarn (https://classic.yarnpkg.com/en/docs/install/#windows-stable)  
Docker: windows(https://www.docker.com/get-started)
ubuntu and based distros `bash curl -fsSL https://get.docker.com -o get-docker.sh sudo sh get-docker.sh `
https://docs.docker.com/engine/install/linux-postinstall/  
Git (https://git-scm.com/downloads)

```bash
git clone https://github.com/pipecraft2/pipecraft
cd pipecraft
yarn run install_pipe
yarn electron:serve
```

