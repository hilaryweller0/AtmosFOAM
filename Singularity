Bootstrap:docker
From:ubuntu:$VERSION

%post
	apt-get update -qq
	apt-get install wget software-properties-common apt-transport-https -y --no-install-recommends
	sh -c "wget -O - http://dl.openfoam.org/gpg.key | apt-key add -"
	add-apt-repository "http://dl.openfoam.org/ubuntu dev" -y
	add-apt-repository "http://atmosfoam-apt.s3-website-eu-west-1.amazonaws.com dev" -y
	apt-get -qq update

	DEBIAN_FRONTEND=noninteractive \
	apt-get install devscripts debhelper openfoam-dev libgdal-dev -y

        DEBIAN_FRONTEND=noninteractive \
	apt-get install atmosfoam-tools -y --allow-unauthenticated --no-install-recommends
	rm -rf /var/lib/apt/lists/*

