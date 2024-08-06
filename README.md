# Shiny app for LAM enrichment analysis of Infinium methylation data

Version 0.1

## Run the app

Requirements: Docker capable computer with minimum 4 CPU threads and 8 GB RAM.

Easiest way to get it working is via the docker image.

```
docker pull mziemann/gmea_app
```

Once you have the image downloaded, you can run the following command:

```
docker run -p 3838:3838 mziemann/gmea_app
```

If you are running this on your local PC, you can access the tool by visiting http://localhost:3838/ in the browser.

If you want to run the app on a server, you will additionally need to configure the firewall to allow port 3838
transfers.

## Report issues

If you are having problems with the app or would like to request a feaure, use the GitHub Issues to raise a ticket.
