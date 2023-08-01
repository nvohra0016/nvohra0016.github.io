---
layout: archive
title: "Publications and Articles"
permalink: /publications/
author_profile: true
---

{% if author.googlescholar %}
  You can also find my articles on <u><a href="{{author.googlescholar}}">my Google Scholar profile</a>.</u>
{% endif %}

{% include base_path %}

Education:
======
* BS, MS in Mathematics, [Indian Institute of Science](http://www.math.iisc.ac.in), Bangalore, India, 2012 - 2018
* MS in Mathematics, [Oregon State University](https://math.oregonstate.edu), Corvallis, 2018 - 2020
* PhD in Mathematics, Oregon State University, 2018 - Present

{% for post in site.publications reversed %}
  {% include archive-single.html %}
{% endfor %}
