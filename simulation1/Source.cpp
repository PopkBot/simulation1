#include<thread>
#include<iostream>
#include <chrono>

#include <SFML/Graphics.hpp>


#include "Rocket.h"





Rocket rocket;
double r = 10;









void thread1() { //процедура потока








	while (true) // если получено сообщ о закрытии окна , то false, в остальное врем€ true , окно ожидает
	{


	}
}


int main() {
	int w = 1080, h = 1080;
	int mouseX = (int)(w / 2);
	int mouseY = (int)(h / 2);
	float mouseSensitivity = 3.0f;
	float speed = 0.1f;
	bool mouseHidden = false;
	bool keyFlag[6] = { false, false, false, false, false, false };
	sf::Vector3f pos = sf::Vector3f(-5.0f, 0.0f, 0.0f);
	sf::Clock clock;
	int framesStill = 1;


	sf::RenderWindow window(sf::VideoMode(h, w), "Tracing", sf::Style::Titlebar | sf::Style::Close); //Titlebar | sf::Style::Close
	window.setFramerateLimit(60);
	window.setPosition(sf::Vector2i(830, 0));

	sf::RenderTexture empityTexture;
	empityTexture.create(h, w);
	sf::Sprite empitySprite = sf::Sprite(empityTexture.getTexture());


	sf::Shader shader;
	shader.loadFromFile("Shader.frag", sf::Shader::Fragment);
	shader.setUniform("u_resolution", sf::Vector2f(h, w));

	/*
	sf::Texture textureGrass;
	if (!textureGrass.loadFromFile("texture/P1000283.jpg")) std::cout Ђ "FOcK U!\n";
	textureGrass.loadFromFile("texture/P1000283.jpg");
	*/

	//std::thread th1(thread1); //объ€вление нового потока
	//th1.detach(); //отсоединение потока от main


	while (window.isOpen()) {
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

		sf::Event event;



		while (window.pollEvent(event))
		{


			if (event.type == sf::Event::MouseWheelScrolled) {
				//std::cout Ђ (double)event.mouseWheelScroll.delta Ђ "\n";
				r += (double)event.mouseWheelScroll.delta;
			}

			if (event.type == sf::Event::Closed)
			{
				window.close();
			}
			if (event.type == sf::Event::MouseMoved)
			{

				if (mouseHidden)
				{
					int mx = event.mouseMove.x - w / 2;
					int my = event.mouseMove.y - h / 2;
					mouseX += mx;
					mouseY += my;
					sf::Mouse::setPosition(sf::Vector2i(w / 2, h / 2), window);
					if (mx != 0 || my != 0) framesStill = 1;
				}
			}
			if (event.type == sf::Event::KeyPressed) {
				if (event.key.code == sf::Keyboard::Escape) {
					window.close();
				}
				if (event.key.code == sf::Keyboard::M) {
					if (!mouseHidden) {
						framesStill = 1;
						window.setMouseCursorVisible(false);
						mouseHidden = true;
					}
					else {
						window.setMouseCursorVisible(true);
						mouseHidden = false;
					}
				}
				/*else if (event.key.code == sf::Keyboard::W) keyFlag[0] = true;
				else if (event.key.code == sf::Keyboard::A) keyFlag[1] = true;
				else if (event.key.code == sf::Keyboard::S) keyFlag[2] = true;
				else if (event.key.code == sf::Keyboard::D) keyFlag[3] = true;
				else if (event.key.code == sf::Keyboard::Space) keyFlag[4] = true;
				else if (event.key.code == sf::Keyboard::LShift) keyFlag[5] = true;*/

				if (event.key.code == sf::Keyboard::Space) { rocket.bThrust = 1; }// rocket.betta = 0;rocket.alpha = 0;rocket.gamma = 0; }
				if (event.key.code == sf::Keyboard::W) { rocket.bTurnY = 1; } //rocket.betta += 0.1f; }
				
				if (event.key.code == sf::Keyboard::S) { rocket.bTurnY = -1; } //rocket.betta -= 0.1f; }
				if (event.key.code == sf::Keyboard::A) { rocket.bTurnZ = -1; } //rocket.gamma += 0.1f; }
				if (event.key.code == sf::Keyboard::D) { rocket.bTurnZ = 1; } //rocket.gamma -= 0.1f;}
				if (event.key.code == sf::Keyboard::Q) { rocket.bTurnX = 1; } //rocket.alpha += 0.1f; }
				if (event.key.code == sf::Keyboard::E) { rocket.bTurnX = -1; } //rocket.alpha -= 0.1f; }
				if (event.key.code == sf::Keyboard::Z) { rocket.wm = sf::Vector3f{0,0,0}; } //rocket.alpha -= 0.1f; }
				//printf("bthrust =%d\n", rocket.bThrust);

			}
			else if (event.type == sf::Event::KeyReleased)
			{
				//if (event.key.code == sf::Keyboard::F) keyFlag[0] = false;
				/*else if (event.key.code == sf::Keyboard::A) keyFlag[1] = false;
				else if (event.key.code == sf::Keyboard::S) keyFlag[2] = false;
				else if

				(event.key.code == sf::Keyboard::D) keyFlag[3] = false;
				else if (event.key.code == sf::Keyboard::Space) keyFlag[4] = false;
				else if (event.key.code == sf::Keyboard::LShift) keyFlag[5] = false;*/
				if (event.key.code == sf::Keyboard::Space)rocket.bThrust = 0;
				if (event.key.code == sf::Keyboard::W)rocket.bTurnY = 0;
				if (event.key.code == sf::Keyboard::S)rocket.bTurnY = 0;
				if (event.key.code == sf::Keyboard::A)rocket.bTurnZ = 0;
				if (event.key.code == sf::Keyboard::D)rocket.bTurnZ = 0;
				if (event.key.code == sf::Keyboard::Q)rocket.bTurnX = 0;
				if (event.key.code == sf::Keyboard::E)rocket.bTurnX = 0;

				if (event.key.code == sf::Keyboard::C) {
					rocket.bAIControl = !rocket.bAIControl;
					if (!rocket.bAIControl)
						{ rocket.bTurnX = 0, rocket.bTurnY = 0, rocket.bTurnZ = 0; }
				}
				
					

			}
		
		}
		rocket.control();
		float mx = ((float)mouseX / w - 0.5f) * mouseSensitivity;
		float my = ((float)mouseY / h - 0.5f) * mouseSensitivity;
		if (mouseHidden) {


			sf::Vector3f dir = sf::Vector3f(0.0f, 0.0f, 0.0f);
			sf::Vector3f dirTemp;
			/*
			if (keyFlag[0]) dir = sf::Vector3f(1.0f, 0.0f, 0.0f);
			else if (keyFlag[2]) dir = sf::Vector3f(-1.0f, 0.0f, 0.0f);
			if (keyFlag[1]) dir += sf::Vector3f(0.0f, -1.0f, 0.0f);
			else if (keyFlag[3]) dir += sf::Vector3f(0.0f, 1.0f, 0.0f);
			dirTemp.z = dir.z * cos(-my) - dir.x * sin(-my);
			dirTemp.x = dir.z * sin(-my) + dir.x * cos(-my);
			dirTemp.y = dir.y;
			dir.x = dirTemp.x * cos(mx) - dirTemp.y * sin(mx) ;
			dir.y = dirTemp.x * sin(mx) + dirTemp.y * cos(mx) ;
			dir.z = dirTemp.z ;
			*/

			dirTemp.z += -r * sin(my);
			dirTemp.x += -r * cos(my) * cos(mx);
			dirTemp.y += -r * cos(my) * sin(mx);
			if (!keyFlag[0]) {
				dir.x = dirTemp.x;
				dir.y = dirTemp.y;
				dir.z = dirTemp.z;
			}
			else {
				dir.x = dirTemp.x * cos(mx) - dirTemp.y * sin(mx);
				dir.y = dirTemp.x * sin(mx) + dirTemp.y * cos(mx);
				dir.z = dirTemp.z;

			}
			pos = dir;
		}
		rocket.dynamic();
		rocket.traektory(shader);
		//rocket.spinning();
		shader.setUniform("u_pos", pos + rocket.rOC);
		shader.setUniform("u_mouse", sf::Vector2f(mx, my));
		shader.setUniform("u_time", clock.getElapsedTime().asSeconds());
		rocket.draw(shader); //shader,textureGrass
		window.draw(empitySprite, &shader);

		window.display();
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		//std::cout Ђ (float)round(round((end - begin).count() /1000)/1000)/1000 Ђ "\n";
		//rocket.dt = (float)round(round((end - begin).count() / 1000) / 1000) / 1000;
		//std::cout Ђ "Time difference = " Ђ std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() Ђ "[µs]" Ђ std::endl;

	}

	return 0;
}