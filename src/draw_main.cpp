#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <util.h> 
#include <graph.h>
#include <kuramoto.h>
#include <time.h>

using namespace std;

#define BOLAS 2

int main()
{
	srand ( time(NULL) );
	Graph<KuramotoOscillator> g;
	KuramotoOscillator::simpleKuramotoNetwork(g, 2, 50);
    // create the window
    sf::RenderWindow window(sf::VideoMode(800, 600), "My window");
	sf::View view(
		sf::Vector2f(0.0, 0.0), 
		sf::Vector2f(200, 150) );
	window.setView(view);
	
    // run the program as long as the window is open
    while (window.isOpen())
    {
		for (int i = 0; i < 2000; ++i) {
        	sf::Event event;
        	while (window.pollEvent(event))
        	{
            	// "close requested" event: we close the window
            	if (event.type == sf::Event::Closed)
                	window.close();
        	}

        	// clear the window with black color
        	window.clear(sf::Color::Black);
			for (int i = 0; i < 2; ++i) {
				sf::CircleShape sp(3);
				sp.setPosition(5*i, g.nodes[i].phase * 30.0);
	//			int color = (57.0 + g.nodes[i].phase)/114.00 * 255;
				if (i) sp.setFillColor( sf::Color::Blue );
				else sp.setFillColor( sf::Color::Red );
				window.draw(sp);
			}
        	// end the current frame
        	window.display();
			g.update();
			sf::sleep( sf::seconds(0.1) );	
		}
    }

    return 0;
}
