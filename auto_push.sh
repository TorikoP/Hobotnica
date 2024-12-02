#!/bin/bash
cd /home/vpal/hobotnica
git add .
git commit -m "Automated backup: $(date)"
git push origin main
