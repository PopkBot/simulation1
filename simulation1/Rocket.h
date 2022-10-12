#pragma once

#include <SFML/Graphics.hpp>
#include<iostream>



//#include"Vec3.h"



#define k 0.2f
#define wind10 0.5f
#define windfactor 0.25f
#define hr 0.68f
#define hC 0.35f


class Rocket
{
public:

	const double kalpha = round(rand() % 360) / 3.1415;

	sf::Vector3f rOC = { 50,1,-100 }, rCB = { hr-hC,0,0 }, rOM, rOA, rCA = {-hC,0,0}, rOB, rAC, rVecCtrl1 = { 0,1,0 }, rVecCtrl2 = { 0,0,1 }, vecThrust, track[100], vecX = { 1,0,0 }, vecY = { 0,1,0 }, vecZ = { 0,0,1 };
	sf::Vector3f startPos = rOC, finishPos;
	float
		ra = 0.15f,
		rb = 0.11f,
		h = hr-hC,
		h2 = 0.5,
		lctrl = rb,
		lvecCtlr = 0.5f;


	int
		trackCount;


	float fire[7] = { 0,0,0,0 ,0,0,0 };
	sf::Vector3f vecCtrl[7];

	sf::Vector3f vecXT;



	float
		dt = 0.01f,
		t = 0,
		m0,
		mf,
		dm,
		m = 4,
		jox = 0.004863f,
		joy = 0.164f,
		joz = 0.164f,
		dax,
		day,
		daz,
		dEox,
		dEoy,
		dEoz,
		alpha = 0,
		betta = 0,
		gamma = 0,
		phi = 0,
		tetta = 0,
		fctrl[7],
		thrust = 100,
		ct = 0.8f,
		cx = 0.2f,
		cy = 0.5f,
		cz = 0.5f,
		sX = 0.01f,
		sY = 0.01f,
		sZ = 0.01f,
		vPre,
		wPre,
		phiPre,
		hIgn=0;

	float yAngle;
	sf::Vector3f yturnVec,
		yturnVecPre;

	sf::Vector3f

		g = {-9.81f,0,0 },// -9.81f
		fvec[7],
		fA,
		vm = { 50,1,0},
		rm,
		wm = { 0,0,0 },
		vA = { 0,0,0 },
		vCD = { h * 0.25f,0,0 }, //våêòîð îò ÖÌ ê ÖÄ
		da = { 0,0,0 };

	__int8
		bThrust = 0,
		bStraifY = 0,
		bStraifZ = 0,
		bTurnX = 0,
		bTurnY = 0,
		bTurnZ = 0;

	bool engRun = false,
		bAIControl = false,
		bhightFound = false,
		bhitGround = false;
	int engIgn = 1;




	void draw(sf::Shader& shader) { //sf::Shader& shader,sf::Texture textureGrass


	
		rOA = rOC + coordTransorm(rCA);
		rOB = rOC + coordTransorm(rCB);




		vecCtrl[1] = coordTransorm(normalize(-cross(vecZ, rCB)) * lctrl ) + rOB;
		vecCtrl[4] = coordTransorm(normalize(-cross(rCB, vecZ)) * lctrl ) + rOB;
		vecCtrl[2] = coordTransorm(normalize(cross(vecY, rCB)) * lctrl) + rOB;
		vecCtrl[3] = coordTransorm(normalize(cross(rCB, vecY)) * lctrl) + rOB;
		vecCtrl[5] = coordTransorm(normalize(cross(vecY, rCB)) * lvecCtlr) + vecCtrl[4];
		vecCtrl[6] = coordTransorm(normalize(cross(rCB, vecY)) * lvecCtlr) + vecCtrl[4];
		vecThrust = coordTransorm(normalize(-vecX * cos(tetta) * cos(phi) + vecY * cos(tetta) * sin(phi) + vecZ * sin(tetta)) * thrust / 20.0f) + rOA;
		
		shader.setUniform("u_thrust", vecThrust);
		shader.setUniform("u_coneROA", rOA);
		shader.setUniform("u_coneROC", rOC);
		shader.setUniform("u_coneROB", rOB);
		shader.setUniform("u_lvecCtlr", lvecCtlr);
		shader.setUniformArray("u_vecCtrl", vecCtrl, 7);
		shader.setUniformArray("u_fire", fire, 7);
		shader.setUniform("u_ra", ra);
		shader.setUniform("u_rb", rb);

		shader.setUniform("u_vecX", coordTransorm(vecX));
		shader.setUniform("u_vecY", coordTransorm(vecY));
		shader.setUniform("u_vecZ", coordTransorm(vecZ));

	

	}

	void spinning() {
		alpha = 0.1f;
		//betta =0.79f;
		//gamma +=0.001f;

		/*if (rOA.x > 20) rOA.x = 0;
		rOA.x += 0.05;
		rOA.y = rOA.x * rOA.x / 20;
		rOA.z = (-rOA.x * 2);
		*/
		/*if (rOA.x * rOA.x + rOA.y * rOA.y + rOA.z * rOA.z < 1000000)
		rOA += coordTransorm(vecFullUnTurn(sf::Vector3f(1, 0,0), alpha,

		betta, gamma));
		else rOA = sf::Vector3f(0, 0, 0);
		*/
		printf("%5f\t%5f\t%5f\t%5f\n", rCB.x, rCB.y, rCB.z, 1);
	//	rAB = qatTurn(rAB, normalize({ 1,1,1 }), alpha);
	//	vecZ = qatTurn(vecZ, normalize({ 1,1,1 }), alpha);
	//	vecY = qatTurn(vecY, normalize({ 1,1,1 }), alpha);
	//	vecX = qatTurn(vecX, normalize({ 1,1,1 }), alpha);
	}


	void dynamic() {

		vecX = normalize(vecX);
		vecY = normalize(vecY);
		vecZ = normalize(vecZ);
		//vm = vecFullTurn(vm, alpha, betta, gamma);

		vA.y = wind10 * sin(kalpha) * pow(abs(rOC.z) / 10, k) + windfactor * round(rand() % 20 - 10) / 10;
		vA.z = wind10 * cos(kalpha) * pow(abs(rOC.z) / 10, k) + windfactor * round(rand() % 20 - 10) / 10;
		vA.x = windfactor * round(rand() % 4 - 2) / 10.0f;


		
		fA.x = cx * 0.635f * abs(dot(vA - vm, vecX)) * dot(vA - vm, vecX) * sX;
		fA.y = cy * 0.635f * abs(dot(vA - vm, vecY)) * dot(vA - vm, vecY) * sY;
		fA.z = cz * 0.635f * abs(dot(vA - vm, vecZ)) * dot(vA - vm, vecZ) * sZ;
		fA = fA.x * vecX + fA.y * vecY + fA.z * vecZ;

		fA = { 0,0,0 };


		/*
		dax = vecFullTurn(g, alpha, betta, gamma).x+(thrust*fire[0]*cos(-tetta)*cos(phi) +fA.x)/m;
		day = vecFullTurn(g, alpha, betta, gamma).y+(thrust*fire[0]*cos(-tetta)*sin(-phi) +fA.y)/m; //+ct*0*(-fire[1]+fire[4])
		daz = vecFullTurn(g, alpha, betta, gamma).z+(thrust*fire[0]*sin(-tetta) +fA.z)/m; //+ct*0*(-fire[2]-fire[5]+fire[3]+fire[6])
		vm.x += dax * dt;
		vm.y += day * dt;
		vm.z += daz * dt;
		*/
		/*
		dax =(thrust * fire[0] * cos(-tetta) * cos(phi)) / m;
		day =(thrust * fire[0] * cos(-tetta) * sin(-phi)) / m; //+ct*0*(-fire[1]+fire[4])
		daz = (thrust * fire[0] * sin(-tetta)) / m; //+ct*0*(-fire[2]-fire[5]+fire[3]+fire[6])
		*/

		sf::Vector3f da1;
		da1 = (normalize(vecX * cos(tetta) * cos(phi) - vecY * cos(tetta) * sin(phi) - vecZ * sin(tetta)) * thrust * fire[0] +fA)/m ;
	
		
		
		
		da = sf::Vector3f{ dot(da1,{1,0,0}), dot(da1,{0,1,0}), dot(da1,{0,0,1}) } + g;
		
		vm += da * dt;

		rOC += coordTransorm(vm * dt);

	
		dEox = (ct * ( -fire[6] + fire[5]) * lctrl ) / jox;//+ fA.y*vCD.z-fA.z*vCD.y	//ct*lctrl/jox;
		wm.x += dEox * dt ;
		
		vecZ = qatTurn(vecZ, normalize(vecX), wm.x*dt);
		vecY = qatTurn(vecY, normalize(vecX), wm.x * dt);
		rCB = qatTurn(rCB, normalize(vecX), wm.x * dt);
		rCA= qatTurn(rCA, normalize(vecX), wm.x * dt);
		if (abs(dot(vecZ, vecY)) > 0.001)
			vecY = cross(vecZ, vecX);
		

		dEoy = (ct * (fire[3]  - fire[2] ) * h - dot(fA, vecZ) * vCD.x + dot(fA, vecX) * vCD.z ) / joy;//+ fA.z*vCD.x-fA.x*vCD.z		ct*h/joy
		wm.y += dEoy * dt ;
		vecZ = qatTurn(vecZ, normalize(vecY), wm.y*dt);
		vecX = qatTurn(vecX, normalize(vecY), wm.y*dt);
		rCB = qatTurn(rCB, normalize(vecY), wm.y*dt);
		rCA = qatTurn(rCA, normalize(vecY), wm.y * dt);
		if (abs(dot(vecZ, vecX)) > 0.001)
			vecX = cross(vecY, vecZ);
		

		dEoz = (ct * (-fire[4] + fire[1]) * h - dot(fA, vecX) * vCD.y + dot(fA, vecY) * vCD.x) / joz;//+ fA.x*vCD.y-fA.y*vCD.x			ct*h/joz
		wm.z += dEoz * dt ;
		vecX = qatTurn(vecX, normalize(vecZ), wm.z*dt);
		vecY = qatTurn(vecY, normalize(vecZ), wm.z*dt);
		rCB = qatTurn(rCB, normalize(vecZ), wm.z*dt);
		rCA = qatTurn(rCA, normalize(vecZ), wm.z * dt);
		if (abs(dot(vecX, vecY)) > 0.001)
			vecY = cross(vecZ, vecX);
		

		alpha += wm.x * dt / 2;
		

	

		


		
		betta -= wm.y * dt/2;
	
		gamma += wm.z * dt/2;
		vecXT = vecXTurned(alpha, betta, gamma);

	

		
		

		t += dt;

	
	}

	

	void wStab() {
		bTurnX = AIturnControl(wm.x, 0, 0, wm.x * wm.x / (2 * ct * lctrl / jox), ct * lctrl / jox, 0.01f, 0.01f); //wm.x stab
		bTurnY = AIturnControl(wm.y, 0, 0, wm.y * wm.y / (2 * ct * h / joy), ct * h / joy, 0.01f, 0.01f); //wm.x stab
		bTurnZ = AIturnControl(wm.z, 0, 0, wm.z * wm.z / (2 * ct * h / joz), ct * h / joz, 0.01f, 0.01f); //wm.x stab
	}


	void control() {


		float wY2=0;
		sf::Vector3f yturnVecContr;
		

		if (!bhightFound) {
			hIgn = findH(vecLength(vm), atanf(vm.x/sqrtf(vm.y * vm.y + vm.z * vm.z))+1.5707963f, abs(g.x), abs(rOC.z), thrust / (m * abs(g.x)), 10000, 0.1f);
			bhightFound = true;
			yAngle = getAngH(vecLength(vm), atanf(vm.x / sqrtf(vm.y * vm.y + vm.z * vm.z)) + 1.5707963f, abs(g.x), abs(rOC.z), hIgn);
			yturnVec =normalize( sf::Vector3f{sqrtf(vm.y*vm.y+vm.z*vm.z)*cos(yAngle)/sin(yAngle),-vm.y,-vm.z});
			float vv0=getVH(vecLength(vm), abs(rOC.z), abs(g.x), hIgn);
			wPre = abs(g.x) * sin(yAngle) / vv0*1.1f;
			phiPre = wPre * wPre / (2 * ct * h / joy) + yAngle;
			vPre = sqrtf(vv0 * vv0 + pow(abs(g.x) * wPre / (ct * h / joy), 2) - 2 * vv0 * cos(yAngle) * abs(g.x) * wPre / (ct * h / joy));


			yturnVecPre = normalize(sf::Vector3f{ sqrtf(vm.y * vm.y + vm.z * vm.z) * cos(phiPre) / sin(phiPre),-vm.y,-vm.z });

			printf("\n\n\nh ignition = %f\tyAng =%5f\twPre = %5f\tphiPre = %5f\tvPre = %5f\n\n\n", hIgn,yAngle,wPre,phiPre,vPre);
			
		}

		float xphi;
		float zphi;
		float yphi;
		float vPhi;
		float vToXAngle;
		sf::Vector3f locN = localVecCoord(normalize(cross(sf::Vector3f{ -1,0,0 }, vm)));
		sf::Vector3f locVm = localVecCoord(normalize(sf::Vector3f{ abs(vm.x),-vm.y,-vm.z }));
		sf::Vector3f locXAxis= localVecCoord(sf::Vector3f{ 1,0,0 });
		sf::Vector3f locyturnVec;
		sf::Vector3f locVmVec= (normalize(sf::Vector3f{ (vm.x),vm.y,vm.z }));

		xphi = acos(dot(normalize(sf::Vector3f{ 0,locN.y,locN.z }), sf::Vector3f{ 0,1,0 }));
		zphi = acos(dot(normalize(sf::Vector3f{ locN.x,locN.y,0 }), sf::Vector3f{ 0,1,0 }));
		yphi = acos(dot(normalize(sf::Vector3f{ locXAxis.x,0,locXAxis.z }), sf::Vector3f{ 1,0,0 }));
		vPhi= acos(dot(normalize(sf::Vector3f{ locVmVec.x,locVmVec.y,locVmVec.z }), sf::Vector3f{ -1,0,0 }));
		vToXAngle= acos(dot(normalize(vm), -vecX));
		
		if ( (vecLength(vm) > vPre & vm.x<0) | (abs(rOC.z)<abs(hIgn) & yphi > yAngle) ) {
			if(dot(vm,vecZ)<0)
				wY2 =- wPre;
			else
				wY2 = wPre;
			locyturnVec = localVecCoord(yturnVec);
			printf("\nF2.1");

		}
		
		else if ( yphi <= yAngle& vm.x < 0) {
			if(dot(vm,vecZ)<0)
				wY2 =- abs(g.x) * sin(vPhi) /(vecLength(vm))*1.2f;
			else
				wY2 = abs(g.x) * sin(vPhi) / (vecLength(vm)) * 1.2f;
			locyturnVec = locVm;
			printf("\nF3");
		}
		else {
			wY2 = 0;
			locyturnVec =localVecCoord( yturnVecPre);
			printf("\nF1");
		}

		printf("\twY2 = % 5f\twy = % 5f\tvToXAngle = %5f", 
			 wY2, wm.y, vToXAngle*180/3.1415f);

		
		

		


	

		
		if(dot(normalize(sf::Vector3f{ 0,locN.y,locN.z }), sf::Vector3f{ 0,0,1 })>=0)
			xphi = acos(dot(normalize(sf::Vector3f{ 0,locN.y,locN.z }), sf::Vector3f{ 0,1,0 }));
		else
			xphi = -acos(dot(normalize(sf::Vector3f{ 0,locN.y,locN.z }), sf::Vector3f{ 0,1,0 }));
		if (abs(xphi) >= 1.57f)
			xphi -= 3.14159265f*abs(xphi)/ xphi;



		if (-dot(normalize(sf::Vector3f{ locN.x,locN.y,0 }), sf::Vector3f{ 1,0,0 }) >= 0)
			zphi = acos(dot(normalize(sf::Vector3f{ locN.x,locN.y,0 }), sf::Vector3f{ 0,1,0 }));
		else
			zphi = -acos(dot(normalize(sf::Vector3f{ locN.x,locN.y,0 }), sf::Vector3f{ 0,1,0 }));
		if (abs(zphi) >= 1.57f)
			zphi -= 3.14159265f * abs(zphi) / zphi;


		/*if (engRun | engIgn == 0) {

			if (-dot(normalize(sf::Vector3f{ locVm.x,0,locVm.z }), sf::Vector3f{ 0,0,1 }) >= 0)
				yphi = acos(dot(normalize(sf::Vector3f{ locVm.x,0,locVm.z }), sf::Vector3f{ 1,0,0 }));
			else
				yphi = -acos(dot(normalize(sf::Vector3f{ locVm.x,0,locVm.z }), sf::Vector3f{ 1,0,0 }));
			if (abs(yphi) >= 1.57f)
				yphi -= 3.14159265f * abs(yphi) / yphi;
		}
		else {
			if (-dot(normalize(sf::Vector3f{ locyturnVec.x,0,locyturnVec.z }), sf::Vector3f{ 0,0,1 }) >= 0)
				yphi = acos(dot(normalize(sf::Vector3f{ locyturnVec.x,0,locyturnVec.z }), sf::Vector3f{ 1,0,0 }));
			else
				yphi = -acos(dot(normalize(sf::Vector3f{ locyturnVec.x,0,locyturnVec.z }), sf::Vector3f{ 1,0,0 }));
			if (abs(yphi) >= 1.57f)
				yphi -= 3.14159265f * abs(yphi) / yphi;
		}
		*/
		if (-dot(normalize(sf::Vector3f{ locyturnVec.x,0,locyturnVec.z }), sf::Vector3f{ 0,0,1 }) >= 0)
			yphi = acos(dot(normalize(sf::Vector3f{ locyturnVec.x,0,locyturnVec.z }), sf::Vector3f{ 1,0,0 }));
		else
			yphi = -acos(dot(normalize(sf::Vector3f{ locyturnVec.x,0,locyturnVec.z }), sf::Vector3f{ 1,0,0 }));
		if (abs(yphi) >= 1.57f)
			yphi -= 3.14159265f * abs(yphi) / yphi;
		


		if (bAIControl ) {
			
			bTurnX = AIturnControl(wm.x, 0, 0, xphi, ct * lctrl / jox, 0.01f,0.01f);
			bTurnZ = AIturnControl(wm.z, 0, 0, zphi, ct * h / joz, 0.01f, 0.01f);
			bTurnY = AIturnControl(wm.y, wY2, 0, yphi, ct * h / joy, 0.001f, 0.01f);


			if ((abs(rOC.z)- hIgn) < 0 & vm.x<0) {
				bThrust = 1;
			
			}

			if (vm.x > 0 & bThrust ==1)
			{
				bThrust = 0;
				printf("\n\n eng stop delta H = %5f", -(rOC.z));
			}

		}

		if (!bhitGround & -rOC.z<0.1f) {
			bhitGround = true;
			finishPos = rOC;
			printf("\n\n\n TOUCH DOWN H = %5f\tV = %5f\tVX = %5f\t AA = %5f\tDistance = %5f", -rOC.z, vecLength(vm), vm.x, getAng(vm.x, sqrtf(vm.y * vm.y + vm.z * vm.z)), vecLength(sf::Vector3f{startPos.x-finishPos.x,startPos.y-finishPos.y,0}));

		}
		
		if (!engRun & bThrust > 0 & engIgn>0) {
			engIgn -= 1;
			engRun = true;
			printf("\n\nEng started H = %5f\tV = %5f\t AA = %5f", -(rOC.z),vecLength(vm),getAng(vm.x,sqrtf(vm.y*vm.y+vm.z*vm.z)));
		}
		if (engRun & bThrust <= 0) {
			engRun = false;
		}
		if (engIgn < 0) engIgn = 0;



		if (engRun) {
			fire[0] = 1;

		


			if (bTurnY == 0 & abs(tetta) > 0) {
				tetta -= abs(tetta) / tetta * 0.02f;
			}
			if (bTurnZ == 0 & abs(phi) > 0) {
				phi -= abs(phi) / phi * 0.02f;
			}
			if (phi > 0.262) phi = 0.262;
			else if (phi < -0.262)phi = -0.262;
			if (tetta > 0.262) tetta = 0.262;
			else if (tetta < -0.262)tetta = -0.262;
		}
		else fire[0] = 0;


		fire[1] = -bTurnZ;
		if (fire[1] < 0) fire[1] = 0;


		fire[2] = bTurnY ;
		if (fire[2] < 0)fire[2] = 0;

		fire[3] = -bTurnY ;
		if (fire[3] < 0)fire[3] = 0;


		fire[4] = bTurnZ;
		if (fire[4] < 0) fire[4] = 0;

		fire[5] = - bTurnX;
		if (fire[5] < 0)fire[5] = 0;

		fire[6] =   bTurnX;
		if (fire[6] < 0)fire[6] = 0;



	}



	void traektory(sf::Shader& shader) {
		float ft = falltime(rOC.z);


		if (ft > 0) {
			trackCount = abs((int)round(ft) / 2) + 20;

			track[0].x = rOC.x - vm.y * ft;
			track[0].y = rOC.y - vm.z * ft;
			track[0].z = 0;
			if (trackCount > 98)trackCount = 98;
			for (int i = 1; i < trackCount; i++) {


				track[i].x = rOC.x - vm.y * ft / (trackCount)*i;
				track[i].y = rOC.y - vm.z * ft / (trackCount)*i;
				track[i].z = rOC.z - vm.x * ft / (trackCount)*i - g.x * powf(ft / (trackCount)*i, 2) / 2;

			}
			shader.setUniformArray("u_track", track, trackCount);
			shader.setUniform("u_trackCount", (float)trackCount);
			


		}

		shader.setUniform("u_Avec", coordTransorm(sf::Vector3f{ dot(fA,{1,0,0}), dot(fA,{0,1,0}), dot(fA,{0,0,1}) })); //vecFullUnTurn(fA, alpha, betta, gamma)	sf::Vector3f{ dot(fA,{1,0,0}), dot(fA,{0,1,0}), dot(fA,{0,0,1}) })
		shader.setUniform("u_CDvec", coordTransorm(normalize(rCB) * vecLength(vCD)));//vecFullUnTurn(vCD,alpha,betta,gamma)
	}


	float falltime(float z0) {
		float ft;
		if (z0 * 0 - rOC.z > 0)
			ft = (vm.x + sqrtf(vm.x * vm.x + 2 * g.x * (-z0 * 0 + rOC.z))) / (-g.x);
		else
			ft = 0;

		return ft;
	}
	sf::Vector3f xUnTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0;
		newVec.y = y0 * cos(alpha0) + z0 * sin(alpha0);
		newVec.z = z0 * cos(alpha0) - y0 * sin(alpha0);
		return newVec;
	}
	sf::Vector3f yUnTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0 * cos(alpha0) - z0 * sin(alpha0);
		newVec.y = y0;
		newVec.z = z0 * cos(alpha0) + x0 * sin(alpha0);
		return newVec;
	}
	sf::Vector3f zUnTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0 * cos(alpha0) + y0 * sin(alpha0);
		newVec.y = y0 * cos(alpha0) - x0 * sin(alpha0);
		newVec.z = z0;
		return newVec;
	}
	sf::Vector3f xTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0;
		newVec.y = y0 * cos(alpha0) - z0 * sin(alpha0);
		newVec.z = z0 * cos(alpha0) + y0 * sin(alpha0);
		return newVec;
	}
	sf::Vector3f yTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0 * cos(alpha0) + z0 * sin(alpha0);
		newVec.y = y0;
		newVec.z = z0 * cos(alpha0) - x0 * sin(alpha0);
		return newVec;
	}
	sf::Vector3f zTurn(sf::Vector3f vec, double alpha0) {
		sf::Vector3f newVec;
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = x0 * cos(alpha0) - y0 * sin(alpha0);
		newVec.y = y0 * cos(alpha0) + x0 * sin(alpha0);
		newVec.z = z0;
		return newVec;
	}
	sf::Vector3f fullUnTurn(double x0, double y0, double z0, double a0, double b0, double g0) {
		sf::Vector3f newVec;
		newVec.x = x0 * cos(b0) * cos(g0) - z0 * sin(b0) + y0 * cos(b0) * sin(g0);
		newVec.y = y0 * (cos(a0) * cos(g0) + sin(a0) * sin(b0) * sin(g0)) - x0 * (cos(a0) * sin(g0) - cos(g0) * sin(a0) * sin(b0)) + z0 * cos(b0) * sin(a0);
		newVec.z = x0 * (sin(a0) * sin(g0) + cos(a0) * cos(g0) * sin(b0)) - y0 * (cos(g0) * sin(a0) - cos(a0) * sin(b0) * sin(g0)) + z0 * cos(a0) * cos(b0);
		return newVec;
	}
	sf::Vector3f vecFullUnTurn(sf::Vector3f vec, double a0, double b0, double g0) {
		sf::Vector3f newVec = vec;

		double x0 = vec.x, y0 = vec.y, z0 = vec.z;

		//newVec.x = z0 * (sin(a0) * sin(g0) + cos(a0) * sin(b0) * cos(g0)) - y0 * (sin(a0) * cos(g0) - cos(a0) * sin(b0) * sin(g0)) + x0 * cos(a0) * cos(b0);
		//newVec.y = y0 * (cos(a0) * cos(g0) + sin(a0) * sin(b0) * sin(g0)) - z0 * (cos(a0) * sin(g0) - cos(g0) * sin(a0) * sin(b0)) + x0 * cos(b0) * sin(a0);
		//newVec.z = y0 * cos(b0) * sin(g0) - x0 * sin(b0) + z0 * cos(b0) * cos(g0);

		newVec = zUnTurn(yUnTurn(xUnTurn(newVec, a0), b0), g0);

		return newVec;
	}
	sf::Vector3f vecFullTurn(sf::Vector3f vec, double a0, double b0, double g0) {
		sf::Vector3f newVec = vec;
		/*
		double x0 = vec.x, y0 = vec.y, z0 = vec.z;
		newVec.x = z0

		* sin(b0) + x0 * cos(b0) * cos(g0) - y0 * cos(b0) * sin(g0);
		newVec.y = x0 * (cos(a0) * sin(g0) + sin(a0) * sin(b0) * cos(g0)) + y0 * (cos(a0) * cos(g0) - sin(a0) * sin(b0) * sin(g0)) - z0 * sin(a0) * cos(b0);
		newVec.z = x0*(sin(a0)*sin(g0)-cos(a0)*sin(b0)*cos(g0))+y0*(sin(a0)*cos(g0)+cos(a0)*sin(b0)*sin(g0))+z0*cos(a0)*cos(b0);
		*/
		//newVec = xTurn(yTurn(zTurn(newVec, g0),b0), a0);
		newVec = xTurn(yTurn(zTurn(newVec, g0), b0), a0);
		return newVec;
	}
	sf::Vector3f coordTransorm(sf::Vector3f vec) {
		sf::Vector3f newVec;
		newVec.x = -vec.y;
		newVec.y = -vec.z;
		newVec.z = -vec.x;
		return newVec;
	}



	sf::Vector3f qatTurn(sf::Vector3f p1, sf::Vector3f u, float a) {


		sf::Vector3f p = p1;

		double t = p.x * u.x + p.y * u.y + p.z * u.z;

		//p.x = p.x * cos(a) + (p.z * u.y - u.z * p.y) * sin(a / 2) * cos(a / 2) + p.x * t * pow(sin(a / 2), 2);
		//p.y = p.y * cos(a) + (p.x * u.z - u.x * p.z) * sin(a / 2) * cos(a / 2) + p.y * t * pow(sin(a / 2), 2);
		//p.z = p.z * cos(a) + (p.y * u.x - u.y * p.x) * sin(a / 2) * cos(a / 2) + p.z * t * pow(sin(a / 2), 2);

		p.x = p1.x * cos(a) + (p1.z * u.y - u.z * p1.y) * sin(a ) + u.x * t * (1-cos(a));
		p.y = p1.y * cos(a) + (p1.x * u.z - u.x * p1.z) * sin(a)  + u.y * t * (1 - cos(a));
		p.z = p1.z * cos(a) + (p1.y * u.x - u.y * p1.x) * sin(a)  + u.z * t * (1 - cos(a));


		p = normalize(p) * vecLength(p1);
		return p;
	}

	sf::Vector3f normalize(sf::Vector3f vec) {

		sf::Vector3f newVec = vec;

		newVec = newVec / vecLength(newVec);
		return newVec;
	}

	sf::Vector3f vecRecoverLength(sf::Vector3f vec, float length) {
		sf::Vector3f newVec = vec;
		newVec = normalize(newVec);
		newVec = newVec * length;
		return newVec;

	}




	float vecLength(sf::Vector3f vec) {
		float length = sqrtf(pow(vec.x, 2) + pow(vec.y, 2) + pow(vec.z, 2));
		//printf("\n\nx = %5f\ty = %5f\tz = %5f\tL = %5f", vec.x, vec.y, vec.z, length);
		return (length);
	}

	sf::Vector3f cross(sf::Vector3f vec1, sf::Vector3f vec2) {
		sf::Vector3f nv1 = vec1, nv2 = vec2, newVec;

		newVec.x = nv1.y * nv2.z - nv1.z * nv2.y;
		newVec.y = nv1.z * nv2.x - nv1.x * nv2.z;
		newVec.z = nv1.x * nv2.y - nv1.y * nv2.x;
		return newVec;
	}

	float dot(sf::Vector3f vec1, sf::Vector3f vec2) {
		float d = vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
		return d;
	}


	sf::Vector3f vecXTurned(float a1, float b1, float g1) {

		sf::Vector3f newVec;
		newVec.x = cos(b1) * cos(g1)-sin(a1)*sin(b1)*sin(g1);
		newVec.y = cos(b1) * sin(g1)+sin(a1)*sin(b1)*cos(g1);
		newVec.z = -cos(a1)*sin(b1);
		normalize(newVec);
		return(newVec);



	}

	void printVecCoord(sf::Vector3f vec) {
		sf::Vector3f newVec;
		float dd = vecX.x * vecY.y * vecZ.z - vecX.x * vecZ.y * vecY.z - vecY.x * vecX.y * vecZ.z + vecZ.x * vecX.y* vecY.z - vecZ.x * vecY.y * vecX.z + vecX.z*vecY.x*vecZ.y;

		//newVec.x = ((vecY.y * vecZ.z - vecZ.y * vecY.z) * vec.x + (vecZ.y * vecX.z - vecX.y * vecZ.z) * vec.y + (vecX.y * vecY.z - vecY.y * vecX.z) * vec.z)/dd;

		//newVec.y = ((vecZ.x * vecY.z - vecY.x * vecZ.z) * vec.x + (vecX.x * vecZ.z - vecZ.x * vecX.z) * vec.y + (vecY.x * vecX.z - vecX.x * vecY.z) * vec.z)/dd;

		//newVec.z = ((vecY.x * vecZ.y - vecZ.x * vecY.y) * vec.x + (vecZ.x * vecX.y - vecX.x * vecZ.y) * vec.y + (vecX.x * vecY.y - vecY.x * vecX.y) * vec.z)/dd;


		newVec.x = ((vecY.y * vecZ.z - vecZ.y * vecY.z) * vec.x + (vecY.z * vecZ.x - vecY.x * vecZ.z) * vec.y + (vecY.x * vecZ.y - vecY.y * vecZ.x) * vec.z) / dd;

		newVec.y = ((vecX.z * vecZ.y - vecX.y * vecZ.z) * vec.x + (vecX.x * vecZ.z - vecX.z * vecZ.x) * vec.y + (vecX.y * vecZ.x - vecX.x * vecY.z) * vec.z) / dd;

		newVec.z = ((vecX.y * vecY.z - vecX.z * vecY.y) * vec.x + (vecX.z * vecY.x - vecX.x * vecY.z) * vec.y + (vecX.x * vecY.y - vecX.y * vecY.x) * vec.z) / dd;


		printf("\n xn = %5f\tyn = %5f\tzn = %5f\tx = %5f\ty = %5f\tz = %5f\ty*z = %5f\n", newVec.x, newVec.y, newVec.z, vec.x, vec.y, vec.z,dot(vecY,vecZ));
		
	}


	sf::Vector3f localVecCoord(sf::Vector3f vec) {
		sf::Vector3f newVec;
		float dd = vecX.x * vecY.y * vecZ.z - vecX.x * vecZ.y * vecY.z - vecY.x * vecX.y * vecZ.z + vecZ.x * vecX.y * vecY.z - vecZ.x * vecY.y * vecX.z + vecX.z * vecY.x * vecZ.y;


		newVec.x = ((vecY.y * vecZ.z - vecZ.y * vecY.z) * vec.x + (vecY.z * vecZ.x - vecY.x * vecZ.z) * vec.y + (vecY.x * vecZ.y - vecY.y * vecZ.x) * vec.z) / dd;

		newVec.y = ((vecX.z * vecZ.y - vecX.y * vecZ.z) * vec.x + (vecX.x * vecZ.z - vecX.z * vecZ.x) * vec.y + (vecX.y * vecZ.x - vecX.x * vecZ.y) * vec.z) / dd;

		newVec.z = ((vecX.y * vecY.z - vecX.z * vecY.y) * vec.x + (vecX.z * vecY.x - vecX.x * vecY.z) * vec.y + (vecX.x * vecY.y - vecX.y * vecY.x) * vec.z) / dd;

		return(newVec);
	}


	
	

	int checkWPhiPos(float w1,float w2,float phi1,float phi2,float EE,float delta,float deltaW) {

		int check=-1;
		//printf("\nw1 = %5f\tw2 = %5f\tph1 = %5f\tph2 = %5f", w1,w2,phi1,phi2);
		if (abs(w1 - w2) < deltaW & abs(phi1 - phi2) < delta) {
			check = 0;
			
		}

		else if (abs(w2) <= deltaW) {
			
			if (
				(w1 >= 0 & (phi1 < phi2 - delta) & (w1 * w1 < w2 * w2 - 2 * EE * (phi1 - phi2 + deltaW))) |
				(w1 <= 0 & (phi1 < phi2 - delta)) |
				(w1 < 0 & (phi1 > phi2 - delta) & (w1 * w1 > w2 * w2 + 2 * EE * (phi1 - phi2 + deltaW))))

				check = 1;
			else if (
				(w1 > 0 & phi1<phi2 + delta & w1 * w1>w2 * w2 - 2 * EE * (phi1 - phi2 - deltaW)) |
				(w1 >= 0 & phi1 > phi2 + delta) |
				(w1 <= 0 & phi1 > phi2 + delta & w1 * w1 < w2 * w2 + 2 * EE * (phi1 - phi2 - deltaW))
				)

				check = 2;
			else if (w1 > w2 + deltaW) check = 3;
			else if (w1 < w2 - deltaW) check = 4;
			else check= - 1;
		}
		else if (w2 > deltaW) {
			
			if (
				(phi1 < phi2 - w2 * w2 / (2 * EE) - delta & (w1 <= 0 | w1 * w1 < w2 * w2 - 2 * EE * (phi1 - phi2 + deltaW))) |
				(w1 > 0 & phi1 > phi2 - w2 * w2 / (2 * EE) - delta & phi1<phi2-delta & (w1 * w1 < w2 * w2 - 2 * EE * (phi1 - phi2 + deltaW) & (w1 * w1 > w2 * w2 + 2 * EE * (phi1 - phi2 + deltaW)))) |
				(w1 <= 0 & phi1 > phi2 - w2 * w2 / (2 * EE) - delta & w1 * w1 > w2 * w2 + 2 * EE * (phi1 - phi2 + deltaW))
				)

				check = 1;

			else if (
				(w1 > 0 & phi1<phi2 + delta & w1 * w1>w2 * w2 - 2 * EE * (phi1 - phi2 - deltaW)) |
				(w1 > 0 & phi1 > phi2 - w2 * w2 / (2 * EE) + delta & phi1 < phi2 + delta & w1 * w1 < w2 * w2 + 2 * EE * (phi1 - phi2 - deltaW)) |
				(w1 <= 0 & phi1 > phi2 - w2 * w2 / (2 * EE) + delta & w1 * w1 < w2 * w2 + 2 * EE * (phi1 - phi2 - deltaW)) |
				(w1 > 0 & phi1 > phi2 + delta)
				)
				check = 2;
			else if (w1 > w2 + deltaW) check = 3;
			else if (w1 < w2 - deltaW) check = 4;
			else check = -1;

		}
		else if (w2<0 & abs(w2)>deltaW) {

			if (
				(w1 < 0 & phi1 < phi2 - delta) |
				(w1 >= 0 & phi1 < phi2 + w2 * w2 / (2 * EE) - delta & w1 * w1 < w2 * w2 - 2 * EE * (phi1 - phi2 + deltaW)) |
				(w1<0 & phi1>phi2 - delta & phi1 < phi2 + w2 * w2 / (2 * EE) - delta & w1 * w1 < w2 * w2 - 2 * EE * (phi1 - phi2 + deltaW)) |
				(w1<w2 & phi1>phi2 - delta & w1 * w1 > w2 * w2 - 2 * EE * (phi1 - phi2 + deltaW))
				)
				check = 1;
			else if (
				(w1 > 0 & phi1<phi2 + w2 * w2 / (2 * EE) + delta & w1 * w1>w2 * w2 - 2 * EE * (phi1 - phi2 - deltaW)) |
				(phi1 > phi2 + w2 * w2 / (2 * EE) + delta & (w1 >= 0 | w1 * w1 < w2 * w2 + 2 * EE * (phi1 - phi2 - deltaW))) |
				(w1<0 & phi1>phi2 + delta & phi1<phi2 + w2 * w2 / (2 * EE) + delta & w1 * w1>w2 * w2 - 2 * EE * (phi1 - phi2 - deltaW) & w1 * w1 < w2 * w2 + 2 * EE * (phi1 - phi2 - deltaW))
				)
				check = 2;
			else if (w1 > w2 + deltaW) check = 3;
			else if (w1 < w2 - deltaW) check = 4;
			else check = -1;

		}
	
		

		return(check);
	}

	int AIturnControl(float w1, float w2, float phi1, float phi2, float EE, float delta,float deltaW) {
		int aiTurn = 0;
		int check = checkWPhiPos(w1, w2, phi1, phi2, EE, delta,deltaW);
		//printf("\n\n check = %d", check);
		switch (check)
		{
		case 0:
			aiTurn = 0;
		//	printf("\n\n\nflag1");
			break;
			
		case 1:case 4:
			aiTurn = -1;
		//	printf("\n\n\nflag2");
			break;
		case 2:case 3:
			aiTurn = 1;
		//	printf("\n\n\nflag3");
			break;
		case -1:
		//	printf("\n===AI FUCKED===\n");
			break;

		default:
		//	printf("\n===AI no output===\n");
			break;
		}
		return(aiTurn);
	}



	float hightCalc(float Rg, float gg, float vv, float aa,int N) {

		

		float cc = (1 + cos(aa)) / (1 - cos(aa));
		float ai;


		float hh=0;
		for (int i = 1; i <= N; i++) {
			ai = aa / N * i;
			hh += pow((48 / (ai * ai * ai * ai - 12 * ai * ai + 48) - 1) * cc, Rg) * pow(sin(aa), 2) * vv * vv / gg * (-72000 * (ai * ai * ai * ai - 12 * ai * ai + 24) / (ai * ai * ai * pow((ai * ai * ai * ai + 20 * ai * ai - 120), 3)));

		}
		hh = hh * aa / N;
		return(hh);




	}

	float hightCalcCOS(float Rg, float gg, float vv, float aa, int N) {




		float cc =pow( (1 + cos(aa)) / (1 - cos(aa)),Rg)* pow(sin(aa), 2) / gg * aa / N;
		
		float ai;


		float hh = 0;
		for (int i = 1; i <= N; i++) {
			ai = aa / N * i;
			hh += pow(((1-cos(ai))/(1+cos(ai))), Rg) * vv * vv  * cos(ai)/pow(sin(ai),3);

		}
		
	
		hh = hh * cc;
		return(hh);




	}






	float getVH(float vv, float h0,float gg,float hh) {

		float vh = sqrtf(vv * vv + 2 * gg * (h0 - hh));
		return (vh);

	}

	float getAngH(float vv, float aa, float gg, float h0, float hh) {

		float bb;
		float vlvh = vv * sin(aa) / (sqrtf(vv * vv * pow(cos(aa), 2) + 2 * gg * (h0 - hh)));

		if (hh == h0 + vv * vv * pow(cos(aa), 2) / (2 * gg))
			bb = 1.570796f;
		else
			bb = atan(vlvh);
		return(bb);

	}


	float getAng(float vx, float vyz) {

		float aang;

		if (vx == 0)aang = 1.57;
		else if (vyz == 0)aang = 0;
		else aang = atan(abs(vyz / vx));

		return(aang);


	}


	float findH(float vv,float aa,float gg,float h0,float Rg,int N,float eps) {

		

		float hmax = h0 + vv * vv * cos(aa)* cos(aa) / (2 * gg);
		printf("\n HMAX = %5f\tvv = %5f\taa = %5f\trg = %5f\th0 = %5f\tgg = %5f", hmax,vv,aa,Rg,h0,gg);
		float hh=hmax/2;
		float hprev = hmax ;
		float hbuf;
		float hcal = hightCalcCOS(Rg, gg, getVH(vv, h0, gg, hh), getAngH(vv, aa, gg, h0, hh), N);
		printf("\nhcal-hh = %5f\thcal = %5f\thh = %5f", (hcal - hh), hcal, hh);


		while (abs(hcal - hh) > eps) {


			/*
			if (hcal < hh) {
				hbuf = hh-abs(hcal-hh)/2;
				hprev = hh;
				hh = hbuf;
				//printf("\n\nf1");
			}
			else
			{
				hbuf = hh + abs(hcal - hh) / 2;
				hprev = hh;
				hh = hbuf;
				//printf("\n\nf2");
			}
			*/
			
			

			hh += (hcal - hh) / 2;
			hcal = hightCalcCOS(Rg, gg, getVH(vv, h0, gg, hh), getAngH(vv, aa, gg, h0, hh), N);
			printf("\nhcal-hh = %5f\thcal = %5f\thh = %5f\tvv = %5f\taa = %5f", (hcal - hh), hcal, hh, getVH(vv, h0, gg, hh), getAngH(vv, aa, gg, h0, hh));
			
			
		}
		printf("\nhcal-hh = %5f", (hcal - hh));

		printf("\n\nV = %5f\taa = %5f", getVH(vv, h0, gg, hh), getAngH(vv, aa, gg, h0, hh));

		//return((hh+hcal)/2);
		return(hh );
	}


	void preTurn(float vv0, float phi0, float EE, float gg) {

		float w0 = gg * sin(phi0) / vv0;
		float phi1 = w0 * w0 / (2 * EE) + phi0;
		float vv1 = sqrtf(vv0 * vv0 + pow(gg * w0 / EE, 2) - 2 * vv0 * cos(phi0) * gg * w0 / EE);







	}






};