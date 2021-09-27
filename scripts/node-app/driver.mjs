

import axios from 'axios'


const config_client = axios.create({baseURL: 'http://127.0.0.1:18086/', timeout: 1000});
const layout_client = axios.create({baseURL: 'http://127.0.0.1:18083/', timeout: 1000});
const shots_client = axios.create({baseURL: 'http://127.0.0.1:18081/', timeout: 1000});

const get_config = async () => {
	const config = await config_client.get('/');
	
	return config.data.config;
};

const get_layout = async () => {
	const layout = await layout_client.post('random/', {});
	return layout.data.layout;
};

const list_shots = async (config, layout) => {
	const shots = await shots_client.post('/', {
		params: {},
		layout
	});
	return shots.data['shot-infos'];
};

const main = async () => {
	const config = await get_config();
	console.log('config', JSON.stringify(config, null, 2));
	const layout = await get_layout();
	console.log('layout', JSON.stringify(layout, null, 2));
	const shots = await list_shots(config, layout);
	console.log('shots', JSON.stringify(shots, null, 2));
}

main();
	
